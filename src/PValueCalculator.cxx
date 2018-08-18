/*==============================================================================
$Id: PValueCalculator.cxx 731612 2016-03-22 15:03:15Z adye $

Class to perform cross-section significance calculation using toys.
This calculator is used by StandardFrequentistDiscovery.
Uses the WorkspaceCalculator base class for control and setup and a
RooStats::HypoTestCalculatorGeneric (eg. RooStats::FrequentistCalculator)
to generate and fit to the toys.

Main parameters:

  fileName         workspace file or list of result files
  wsName           workspace name
  modelSBName      ModelConfig object name
  dataName         observed dataset

  calculatorType = 0 Frequentist calculator
                 = 1 Hybrid calculator
                 = 2 Asymptotic calculator

  testStatType   = 0 Simple Likelihood Ratio (LEP)
                 = 1 Ratio of Profile Likelihood (Tevatron)
                 = 2 Profile Likelihood Ratio two-sided
                 = 3 Profile Likelihood Ratio one-sided           (i.e.   0  if mu_hat >= mu)
                 = 4 Maximum Likelihood Estimate (mu_hat)
                 = 5 Profile Likelihood Ratio one-sided discovery (i.e.   0  if mu_hat <= mu)
                 = 6 Profile Likelihood Ratio signed discovery    (i.e. -q0  if mu_hat <  mu)
                 = 7 Profile Likelihood Ratio signed              (i.e. -qmu if mu_hat >  mu)

  ntoys:           number of toys to use

  plotResult:      plot result of test (TS distribution)
  writeResult:     write result (default is true)

Author: Tim Adye, based on $ROOTSYS/tutorials/roostats/StandardFrequentistDiscovery.C from ROOT 5.34.00.

==============================================================================*/

#include "Utilities/WorkspaceCalculatorConfig.h"
#include "PValueCalculator.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TList.h"
#include "TKey.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TText.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStopwatch.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"

#include "RooStats/HypoTestCalculatorGeneric.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"

#include "RooStats/HypoTestResult.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::min;
using std::max;
using std::isnan;

void PValueCalculator::Init() {
   bPOIdefault       = 0.0;
  sbPOIdefault       = 1.0;
  resultFileName     = "p0.root";
  resultName         = "HypoTestCalculator_result";
  fitLeadbetter      = false;  // fit Leadbetter term to q0 distribution
  runToys            = true;

  res                = 0;
}


void PValueCalculator::SetDefaults() {
  if (testStatType<0) testStatType= 2; // test statistic type (2=two-sided)
  // Do data TS for the first job.
  // Also if using importance sampling, so PLTS's NLL is set up before first toy.
  if (jobNum<=0 || samplerType>=0) optimize &= (~kSkipDataTS);

  WorkspaceCalculator::SetDefaults();
}


void PValueCalculator::ShowParms() const
{
     cout << GetName()
          << " wsName=\""         << wsName      << "\""
          << ", modelSBName=\""   << modelSBName << "\""
          << ", modelBName=\""    << modelBName  << "\""
          << ", dataName=\""      << dataName    << "\"";
     PrintOptimize (optimize & ~(kAdjustRanges|kFixStatError|kMinosData), ", ");
     if (optimize & kAdjustRanges)   cout << ", AdjustRanges "<<newParameterRanges<<" sigma";
     if ((optimize & kFixStatError) && fixStatErrorMin>=0)
                                     cout << ", FixStatError "<<100.0*fixStatErrorMin<<"%";
     if (optimize & kMinosData) {
                                     cout << ", MinosData";
       if (minosSetSize>0)           cout << " (set " << minosSetNum << "/" << minosSetSize << ")";
     }
     if (detailedOutput)             cout << ", DetailedOutput";
     if (dropNegative)               cout << ", dropNegative";
     cout << ", ntoys="           << ntoys;
     if (skipToys>0)
     cout << ", skipToys="        << skipToys;
     if (ntoysAlt>0)
     cout << ", ntoysAlt="        << ntoysAlt;
     if (nworkers>0)
     cout << ", nworkers="        << nworkers;
     if (invMass>0.0)
     cout << ", invMass="         << invMass;
     if (!poimin.empty())
     cout << ", poimin="          << Join(poimin);
     if (!poimax.empty())
     cout << ", poimax="          << Join(poimax);
     cout << endl;
}


int PValueCalculator::RunOverWorkspace()
{
  // run the calculator
  int     ok= SetupWorkspace();
  if (ok) ok= SetupMinimizer();
  if (ok) ok= SetupModel();
  if (ok) ok= SetupInitialParms();
  if (ok) ok= SetupTestStat();
  if (ok) ok= SetupSampler();
  if (ok) ok= SetupCalculator();
  if (ok) ok= RunCalculatorProof();
  if (ok) ok= GetTestStatInfo();
#ifdef USE_ToyMCImportanceSampler
  if (samplerType>=0) AddWeights (statsTree, res);
#endif
  DeleteCalculator();   // don't use from here on
  return ok;
}

int PValueCalculator::RunCalculator()
{
   // generate toys, sampling distribution, and get result
   TStopwatch tw;
   res = calc->GetHypoTest();
   if (verbose>=0) {
     cout << "Time to generate toys: "; tw.Print();
   }
   if (!res) return 0;
   if (RooStats::SamplingDistribution* samp= res->GetNullDistribution()) samp->SetTitle("null / b-only");
   if (RooStats::SamplingDistribution* samp= res->GetAltDistribution())  samp->SetTitle("alt / s+b");
   return 1;
}


void PValueCalculator::AdjustResults (bool resultWasRead)
{
  if (!res) return;
  WorkspaceCalculator::AdjustResults (resultWasRead);
  saved.Add (res);
}


//==============================================================================
// Other methods
//==============================================================================

TTree* PValueCalculator::GetLimits()
{
  if (!res) return 0;

  Int_t nToys=0, nToysSig=0;
  RooStats::SamplingDistribution* null= res->GetNullDistribution();
  if (null) nToys= null->GetSize();
  RooStats::SamplingDistribution* alt=  res->GetAltDistribution();
  if (alt) nToysSig= alt->GetSize();

  Double_t tsdata     = res->GetTestStatisticData();
  Double_t pvalue     = res->NullPValue();
  Double_t pvalue_err = res->NullPValueError();
  Double_t exp_err=0.0, q0exp=NaN, exp = ExpectedPValue (res, &exp_err, &q0exp);
  int nFail= 0;
  if (statsTree && testStatType==6)
    nFail= statsTree->GetEntries (Form("!is_data && ml_cond-ml<(%g)", negVal));
  else
    nFail= NumFailures(res->GetNullDistribution());

  Double_t pverr= 0.0, exp_err2=0.0;
  if (doneRead && !haveWeights) pverr= ResampleErrors (res, &exp_err2, nToysResample);

  int npoi= bPOI.size();
  bool twoSided= asimovSig>0.0 || (testStatType==2 && !(npoi==1 && poimin[0]>=bPOI[0])); // really two-sided
  double sides= twoSided ? 2.0 : 1.0;

  if (verbose>=0)
    cout << nToys << " null toys, " << nFail << " fit failures" << endl;
  cout << "Null p-value " << pvalue << " +/- " << pvalue_err;
  if (pverr>0.0) cout << " +/- " << pverr;
  if (pvalue>0.0) {
    double sig= -ROOT::Math::gaussian_quantile(pvalue/sides,1.0);
    cout   << Form(", significance %.3f +/- %.3f",sig,
                   fabs (ROOT::Math::gaussian_quantile((pvalue-pvalue_err)/sides,1.0) + sig));  // larger error for p<0.5
    if (pverr>0.0)
      cout << Form(" +/- %.3f",
                   fabs (ROOT::Math::gaussian_quantile((pvalue-pverr     )/sides,1.0) + sig));
  }
  cout << endl;

  double q= 2.0*tsdata, p, Z;
  if (npoi==1 && asimovSig>0.0) {
    TF1* chi2= AsymTF1 (testStatType, 1);
    p= chi2->Eval(q);
    delete chi2;
    Z= -ROOT::Math::gaussian_quantile (0.5*p, 1.0);
  } else if (npoi==1) {  // shortcut
    Z= sqrt (fabs(q));
    if (q<0.0) Z= -Z;    // if twoSided, gives p>1 for "impossible" q,Z<0 (unlike below, which gives p=1)
    p= sides * ROOT::Math::gaussian_cdf(-Z);
  } else if (twoSided) {
    p= ROOT::Math::chisquared_cdf_c (q, double(npoi));
    Z= -ROOT::Math::gaussian_quantile (0.5*p, 1.0);
  } else {
    p= 0.5*ROOT::Math::chisquared_cdf_c (fabs(q), double(npoi));
    if (q<0.0) p= 1.0-p;
    Z= -ROOT::Math::gaussian_quantile (p, 1.0);
  }

  cout << "Asymptotic p-value " << p << ", significance " << Z << " ("<<npoi<<"DoF, "<<sides<<"-sided";
  if (npoi==1 && asimovSig>0.0) cout << ", Asimov sigma "<<asimovSig;
  cout << ")" << endl;
  if (exp>=0.0) {
    cout << "Expected p-value " << exp << " +/- " << exp_err;
    if (exp_err2>0.0) cout << " +/- " << exp_err2;
    if (!isnan(q0exp)) cout << " - median q"<<bPOI[0]<<"|"<<sbPOI[0]<<" " << q0exp << endl;
    cout << endl;
  }

  TF1* quant= AsymTF1 (testStatType, 2);
  cout << "Likelihood limits: ";
  for (int sig=1; sig<=3; sig++) {
    if (sig>1) cout << ", ";
    double pv= ROOT::Math::gaussian_cdf(-sig);
    cout << Form("%.2f%%: %.3f", 100.0*(1.0-2.0*pv), quant->Eval(sides*pv));
  }
  cout << endl;
  delete quant;

  UInt_t initialSeed= seed;
  TTree* bandTree= new TTree ("band", GetName());
  bandTree->SetDirectory(0);   // memory resident ntuple
  bandTree->Branch ("obsUL",        &pvalue);
  bandTree->Branch ("tsdata",       &tsdata);
  bandTree->Branch ("obsUL_err",    &pvalue_err);
  bandTree->Branch ("invMass",      &invMass);
  bandTree->Branch ("nToys",        &nToys);
  bandTree->Branch ("nToysSig",     &nToysSig);
  bandTree->Branch ("nFail",        &nFail);
  bandTree->Branch ("muhat",        &poihat);
  bandTree->Branch ("muhat_err",    &poierr);
  bandTree->Branch ("muMin",        &poimin);
  bandTree->Branch ("muMax",        &poimax);
  bandTree->Branch ("muBkg",        &bPOI);
  bandTree->Branch ("muSig",        &sbPOI);
  bandTree->Branch ("tstype",       &testStatType);
  bandTree->Branch ("wsfile",       &wsfile);
  bandTree->Branch ("bandMedian",   &exp);
  bandTree->Branch ("bandMedian_err",&exp_err2);
  bandTree->Branch ("optimize",     &optimize);
  bandTree->Branch ("seed",         &initialSeed);
  bandTree->Branch ("args",         &cmdArgs);
  bandTree->Branch ("poiName",      &poiName);
  bandTree->Fill();   // just one entry so we can combine later
  bandTree->ResetBranchAddresses();

  saved.Add (bandTree);
  return bandTree;
}


int PValueCalculator::ReadResult (TFile* f, const char* objname)
{
  RooStats::HypoTestResult* radd = 0;
  f->GetObject (objname, radd);
  if (!radd) return 1;
  if (verbose>=1) PrintResult (radd);
  if (!res) res= dynamic_cast<RooStats::HypoTestResult*>(radd->Clone());
  else      res->Append(radd);
  delete radd;
  return 0;
}


int PValueCalculator::ReadTree (TFile* f, UInt_t& thisSeed)
{
  thisSeed= 0;
  TTree* bandTree = 0;
  f->GetObject ("band", bandTree);
  if (!bandTree) return 0;
  Int_t nent = bandTree->GetEntries();
  if (nent != 1) {
    cerr << "results file contains " << nent << " entries" << endl;
    if (!force) return 1;
    if (nent==0) return 0;
  }

  vector<Double_t> poimin1, poimax1, bPOI1, sbPOI1;
  vector<string>* poiName1=0;
  Double_t invMass1=0.0, poihat1=0.0, poierr1=0.0;
  Int_t testStatType1=-1;
  string* wsfile1=0;
  ULong64_t optimize1=0;
  bandTree->SetBranchAddress ("invMass", &invMass1);
  if (bandTree->GetBranch("muhat"))     bandTree->SetBranchAddress ("muhat",     &poihat1);
  if (bandTree->GetBranch("muhat_err")) bandTree->SetBranchAddress ("muhat_err", &poierr1);
  if (bandTree->GetBranch("tstype"))    bandTree->SetBranchAddress ("tstype",    &testStatType1);
  if (bandTree->GetBranch("wsfile"))    bandTree->SetBranchAddress ("wsfile",    &wsfile1);
  SetBranchAddress (bandTree, "muMin", &poimin1, true, NaN);
  SetBranchAddress (bandTree, "muMax", &poimax1, true, NaN);
  SetBranchAddress (bandTree, "muBkg", &bPOI1,   true, NaN);
  SetBranchAddress (bandTree, "muSig", &sbPOI1,  true, NaN);
  if (bandTree->GetBranch("optimize"))  bandTree->SetBranchAddress ("optimize",  &optimize1);
  if (bandTree->GetBranch("seed"))      bandTree->SetBranchAddress ("seed",      &thisSeed);
  if (bandTree->GetBranch("poiName"))   bandTree->SetBranchAddress ("poiName",   &poiName1);
  bool haveVec= (poimin1.size()==1);
  bandTree->GetEntry(0);
  if (wsfile1 && !wsfile1->empty()) {
    if       (wsfile.empty())
      wsfile= *wsfile1;
    else if (*wsfile1!=wsfile) {
      cerr << "results from different workspace files: " << wsfile << " and " << *wsfile1 << endl;
      if (!force) return 2;
    }
  }
  if (invMass1>0.0) {
    if      (invMass<=0.0)
      invMass= invMass1;
    else if (invMass1!=invMass) {
      cerr << "results at different masses: " << invMass << " and " << invMass1 << " GeV" << endl;
      if (!force) return 3;
    }
  }
  if (testStatType1>=0) {
    if (testStatType<0)
      testStatType= testStatType1;
    else if (testStatType1!=testStatType) {
      cerr << "results at different test statistic types: " << testStatType << " and " << testStatType1 << endl;
      overrideTestStatType= testStatType1;
      if (!force) return 4;
    }
  }
  if  (!poimin1.empty() && !isnan(poimin1[0]) &&
       !poimax1.empty() && !isnan(poimax1[0])) {
    if (poimin.empty()) poimin= poimin1;
    if (poimax.empty()) poimax= poimax1;
    if ((poimin1!=poimin && !(haveVec && poimin.size()>1 && poimin1.size()==1 && poimin1[0]==poimin[0])) ||
        (poimax1!=poimax && !(haveVec && poimax.size()>1 && poimax1.size()==1 && poimax1[0]==poimax[0]))) {
      cerr << "results with different POI ranges: " << Join(poimin,poimax)<< " and " << Join(poimin1,poimax1) << endl;
      if (!force) return 5;
    }
  }
  if (! bPOI1.empty() && !isnan( bPOI1[0])) {
    if ( bPOI.empty())  bPOI=  bPOI1;
    if ( bPOI1!= bPOI && !(haveVec &&  bPOI.size()>1 &&  bPOI1.size()==1 &&  bPOI1[0]== bPOI[0])) {
      cerr << "results with different Bkg POI settings: " << Join( bPOI,":")<< " and " << Join( bPOI1,":") << endl;
      if (!force) return 6;
    }
  }
  if (!sbPOI1.empty() && !isnan(sbPOI1[0])) {
    if (sbPOI.empty()) sbPOI= sbPOI1;
    if (sbPOI1!=sbPOI && !(haveVec && sbPOI.size()>1 && sbPOI1.size()==1 && sbPOI1[0]==sbPOI[0])) {
      cerr << "results with different S+B POI settings: " << Join(sbPOI,":")<< " and " << Join(sbPOI1,":") << endl;
      if (!force) return 6;
    }
  }
  if (poierr==0.0 && poierr1!=0.0) {
    poihat= poihat1;
    poierr= poierr1;
  }
  optimize |= optimize1;
  if (poiName1 && !poiName1->empty()) {
    if (poiName.empty()) poiName= *poiName1;
    if (*poiName1!=poiName) {
      cerr << "results with different POIs: " << Join(poiName,":")<< " and " << Join(*poiName1,":") << endl;
      if (!force) return 6;
    }
  }
  bandTree->ResetBranchAddresses();
  delete bandTree;
  return 0;
}


Double_t PValueCalculator::ExpectedPValue (RooStats::HypoTestResult* r, Double_t* experr, Double_t* q0exp, double nsig)
{
  RooStats::SamplingDistribution* alt= r->GetAltDistribution();
  if (!(alt  && alt ->GetSize())) return -1.0;
  RooStats::SamplingDistribution* null= r->GetNullDistribution();
  if (!(null && null->GetSize())) return -1.0;
  double pv= (nsig==0.0) ? 0.5 : ROOT::Math::normal_cdf(nsig,1);
  Double_t q0exp1= alt->InverseCDFInterpolate(pv);
  if (q0exp) *q0exp= q0exp1;
  if (experr) {
    *experr=0.0;
    return null->IntegralAndError (*experr, q0exp1, RooNumber::infinity(), kTRUE, kTRUE, kTRUE);
  } else {
    return null->Integral                  (q0exp1, RooNumber::infinity(), kTRUE, kTRUE, kTRUE);
  }
}


Double_t
PValueCalculator::ResampleErrors (const RooStats::HypoTestResult* rold, Double_t* experr, int nToys)
{
  Double_t obs=0.0, obs2=0.0, exp=0.0, exp2=0.0;
  if (nToys<2) return 0.0;
  ::Info("PValueCalculator::ResampleErrors","resample p-value distributions %d times to estimate errors", nToys);
  for (int i=0; i<nToys; i++) {
    RooStats::HypoTestResult* r= ResampleToy(rold);
    Double_t v= r->NullPValue();
    obs  += v;
    obs2 += v*v;
    if (experr) {
      Double_t e= ExpectedPValue (r);
      exp  += e;
      exp2 += e*e;
    }
    delete r;
  }
  if (experr) *experr= sqrt (fabs (exp2 - (exp*exp)/nToys) / (nToys-1));
  return               sqrt (fabs (obs2 - (obs*obs)/nToys) / (nToys-1));
}


RooStats::HypoTestResult*
PValueCalculator::ResampleToy (const RooStats::HypoTestResult* rold)
{
  RooStats::HypoTestResult* r= dynamic_cast<RooStats::HypoTestResult*>(rold->Clone());
  RooStats::SamplingDistribution* null= ResampleToy (r->GetNullDistribution());
  if (null) r->SetNullDistribution(null);
  RooStats::SamplingDistribution* alt=  ResampleToy (r->GetAltDistribution());
  if (alt)  r->SetAltDistribution(alt);
  return r;
}


void
PValueCalculator::PrintResult (RooStats::HypoTestResult* r)
{
  r->Print();
  cout << "Result";
  if (!wsfile.empty()) cout << " for workspace file " << wsfile;
  if (invMass>0.0) cout << " at " << invMass << " GeV";
  if (!poiName.empty()) cout << " " << Join(poiName);
  else                  cout << " mu";
  if (poierr!=0.0) cout << " fitted = " << poihat << " +/- " << poierr;
  if (!bPOI.empty() || !sbPOI.empty()) cout << " null="<<Join(bPOI)<<" alt="<<Join(sbPOI);
  if (!poimin.empty() || !poimax.empty()) cout << " in range "<<Join(poimin,poimax);
  if (testStatType>=0) cout << " and test stat #" << testStatType;
  cout <<endl;
  PrintSamplingDistribution (r->GetNullDistribution(), "null");
  PrintSamplingDistribution (r->GetAltDistribution(),  "alt ");
}


int PValueCalculator::LimitSamples()
{
  RooStats::SamplingDistribution* null= LimitSamples (res->GetNullDistribution(), ntoys,    jobNum);
  if (null) res->SetNullDistribution(null);
  RooStats::SamplingDistribution* alt=  LimitSamples (res->GetAltDistribution(),  ntoysAlt, jobNum);
  if (alt)  res->SetAltDistribution(alt);
  return (null || alt) ? 1 : 0;
}


int PValueCalculator::DropBadSamples()
{
  haveWeights= 0;
  RooStats::SamplingDistribution* null= DropBadSamples (res->GetNullDistribution());
  if (null) res->SetNullDistribution(null);
  RooStats::SamplingDistribution* alt=  DropBadSamples (res->GetAltDistribution());
  if (alt)  res->SetAltDistribution(alt);
  return (null || alt) ? 1 : 0;
}

int PValueCalculator::ApplyCuts()
{
  RooStats::SamplingDistribution* null= ApplyCuts (res->GetNullDistribution());
  if (null) res->SetNullDistribution(null);
  RooStats::SamplingDistribution* alt=  ApplyCuts (res->GetAltDistribution(), true);
  if (alt)  res->SetAltDistribution(alt);
  return (null || alt) ? 1 : 0;
}




void PValueCalculator::PlotTSdiff (TPad* canvas, TH1D* he, int tsType, Double_t obsTS, const TF1* fit)
{
  if (tsType==2 && bPOI.size()==1 && poimin[0]>=bPOI[0]) tsType=5;   // like single-sided
  int nbin=he->GetNbinsX();
  double xlo=he->GetXaxis()->GetXmin(), xhi=he->GetXaxis()->GetXmax();
  TString tsvar=tsName(tsType), chi2label=chi2Name(tsType);

  TH1D* hd= new TH1D(Form("%sdiff",he->GetName()), Form ("%s;%s;Entries / %s", he->GetTitle(), tsvar.Data(), chi2label.Data()), nbin, xlo, xhi);
  hd->SetStats(kFALSE);
  hd->SetLineWidth(2);
  hd->SetLineColor(he->GetLineColor());

  Double_t maxd=4.0, ne= he->GetSum();
  if (ne>0.0) {
    TF1* pvfun= AsymTF1 (tsType, 1, he);
    for (int i=1; i<=nbin; i++) {
      double e=  he->GetBinContent(i);
      double ee= he->GetBinError(i);
      if (e==0.0) continue;
      double x0= hd->GetBinLowEdge(i), x1=hd->GetBinLowEdge(i+1);
      Double_t f= pvfun->Eval(x0) - pvfun->Eval(x1);
      if (f==0.0) continue;
      double d= e/(ne*f);
      if (verbose>=2)
        cout << Form("%g:%g: f=%g e=%g e/f=%g +/- %g",x0,x1,f,e,d,d/(sqrt(ne)*ee)) << endl;
      double de= d/ee;
      hd->SetBinContent (i, d);
      hd->SetBinError   (i, de);
      // Adjust histogram maximum so as not to miss significant points with Draw opt "e0"
      if (d-(2.0*de)>=1.0 && 1.1*d>maxd) maxd=1.1*d;
    }
    delete pvfun;
  }
  hd->SetMinimum(0.0);
  hd->SetMaximum(maxd);
  TLine* zline= new TLine (xlo, 1.0, xhi, 1.0);
  zline->SetLineColor(kRed);
  zline->SetLineWidth(2);

  TLine* line= 0;
  if (!isnan(obsTS) && obsTS<=xhi) {
    line= new TLine (obsTS, 0.0, obsTS, maxd);
    line->SetLineWidth(3);
  }

  TF1* fitd= 0;
  if (fit) {
    fitd= new TF1 ("fitd", "2.0*[0]*[3]*((1.0-[1]) + [4]*[1]*[2]*exp([2]*([5]-x))/ROOT::Math::chisquared_pdf(x,1.0))", fit->GetXmin(), xhi);
    fitd->SetParameters (fit->GetParameters());
    fitd->SetNpx(10*nbin);
    fitd->SetLineColor(kMagenta);
  }

  hd->Draw("e0");
  zline->Draw();
  if (line) line->Draw();
  if (fitd) fitd->Draw("lsame");

  PrintPlot (canvas);

  delete fitd;
  delete hd;
  delete line;
  delete zline;
}

void PValueCalculator::PlotPValue (TPad* canvas, TH1* h, RooStats::HypoTestResult* r, int tsType, const TF1* fit)
{
  if (tsType==2 && bPOI.size()==1 && poimin[0]>=bPOI[0]) tsType=5;   // like single-sided
  int nbin=h->GetNbinsX();
  double xhi=h->GetXaxis()->GetXmax();
  TString tsvar=tsName(tsType), chi2label=chi2Name(tsType);

  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle(tsvar);
  h->GetYaxis()->SetTitle("p-value");

  Double_t obsTS = r->GetTestStatisticData(), q0= 2.0*obsTS;
  for (Int_t i= 0; i<=nbin; i++) {
    Double_t x= h->GetBinCenter(i);
    r->SetTestStatisticData(0.5*x);
    h->SetBinContent (i, r->NullPValue());
    h->SetBinError   (i, r->NullPValueError());
  }
  if (verbose>=0) {
    for (Double_t x= 0.0; x<xhi; x += 1.0) {
      if (x==0.0) x= 1e-6;
      r->SetTestStatisticData(0.5*x);
      cout << Form("With q0,obs %.1f p-value = %.3g +/- %.3g, Z = %.2f (asymptotic %.2f)",
                   x, r->NullPValue(), r->NullPValueError(),
                   -ROOT::Math::gaussian_quantile(r->NullPValue(),1.0),
                   (x<0?-1:1)*sqrt(fabs(x))) << endl;
      if (r->NullPValue()<=0.0) break;
    }
  }
  r->SetTestStatisticData(obsTS);

  TF1* cdf= AsymTF1 (tsType, 1);
  Double_t obsp= cdf->Eval(xhi);
  if (obsp>0.0) h->SetMinimum(obsp);

  TF1* fitcdf= 0;
  if (fit) {
    fitcdf= new TF1 ("fitcdf", "[0]*[3]*((1.0-[1])*2.0*ROOT::Math::gaussian_cdf(-sqrt(x)) + [4]*[1]*exp([2]*([5]-x)))", fit->GetXmin(), xhi);
    fitcdf->SetParameters (fit->GetParameters());
    fitcdf->SetNpx(10*nbin);
    fitcdf->SetLineColor(kMagenta);
    if (!isnan(q0) && verbose>=0) {
      Double_t obspp= fitcdf->Eval(q0);
      cout << "Fitted p-value at q0 "<<q0<<": p0 = "<<obspp<<", Z = "<<-ROOT::Math::gaussian_quantile(obspp,1.0)<<endl;
    }
  }

  TLine* line= 0;
  if (!isnan(q0) && q0<=xhi) {
    Double_t y= h->GetBinContent (h->FindBin(q0));
    if (y<0.01*h->GetMaximum()) y= 0.1*h->GetMaximum();
    line= new TLine (q0, 0.0, q0, y);
    line->SetLineWidth(3);
  }

  Double_t xoff=0.99-canvas->GetRightMargin(), yoff=0.99-canvas->GetTopMargin();
  TLegend* leg= new TLegend (xoff-0.18, yoff-0.12, xoff, yoff);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
            leg->AddEntry (h,      "Toys",                                "L");
  if (line) leg->AddEntry (line,   "Observed",                            "L");
            leg->AddEntry (cdf,    chi2label,                             "L");
  if (fit)  leg->AddEntry (fitcdf, Form("%s + exp fit",chi2label.Data()), "L");

  h->Draw("e0");
  cdf->Draw("lsame");
  if (fitcdf) fitcdf->Draw("lsame");
  h->Draw("e0same");
  if (line) line->Draw();
  leg->Draw();

  PrintPlot (canvas);

  delete line;
  delete fitcdf;
}

Double_t PValueCalculator::GetVal (TTree* tree, const char* name, Long64_t entry, Double_t defval)
{
  if (entry>=0) tree->GetEntry(entry);
  TLeaf* leaf= tree->FindLeaf(name);
  if (!leaf) return defval;
  return leaf->GetValue();
}

int PValueCalculator::TopCorr (const std::vector<string>& npnames,
                                  const int npoi,
                                  std::vector<int>& itop,
                                  std::vector<int>& jtop,
                                  std::vector<Double_t>& topcor,
                                  TH2D** hcor_ret,
                                  const Double_t maxCorrPOI,
                                  const Double_t maxCorr,
                                  const int maxTop)
{
  // Find top pairs of correlated nuisance parameters
  int nnp= npnames.size();
  if (nnp<2)  return 0;
  const Long64_t nent= statsTree->GetEntries();
  if (nent<3) return 0;

  Int_t is_data=0, is_alt=0;
  Double_t weight= 1.0, sumw=0.0;
  TVectorD val(nnp), sumi(nnp);
  TMatrixD sumij(nnp,nnp), cov(nnp,nnp), cor(nnp,nnp);
  Long64_t nsel=0;
  TTreeFormula* select=0;
  if (cut.Length()>0) {
    select = new TTreeFormula ("selection", cut, statsTree);
    if (!select->GetNdim()) {
      cerr << "bad cut: " << cut << endl;
      delete select; select=0;
    }
  }
  statsTree->SetBranchStatus("*",0);
  SetBranchAddress (statsTree, "is_data", &is_data);
  SetBranchAddress (statsTree, "is_alt",  &is_alt, true);
  SetBranchAddress (statsTree, "weight",  &weight, true);
  for (int i=0; i<nnp; i++)
    SetBranchAddress (statsTree, npnames[i].c_str(), &val[i]);
  if (select)
    for (Int_t i=0, n=select->GetNcodes(); i<n; i++)
      select->GetLeaf(i)->GetBranch()->SetStatus(1);
  for (Long64_t entry=0; entry<nent; entry++) {
    statsTree->GetEntry(entry);
    if (is_data)     continue;
    if (is_alt)      continue;
    if (weight==0.0) continue;
    if (plotMuHatPositive && val[0]<0.0) continue;
    if (select) {
      Int_t i=0, ndata=select->GetNdata();
      for (; i<ndata; i++) {
        if (select->EvalInstance(i) != 0.0) break;
      }
      if (i>=ndata) continue;
    }
    nsel++;
    sumw += weight;
    for (int i=0; i<nnp; i++) {
      Double_t ival= weight*val[i];
      sumi[i] += ival;
      for (int j=0; j<=i; j++)
        sumij(i,j) += weight*ival*val[j];
    }
  }
  statsTree->ResetBranchAddresses();
  statsTree->SetBranchStatus("*",1);
  delete select;
  if (nsel<3) return 0;

  // Could multiply cov(i,j) by nsel/(nsel-1) for sample covariances, but this cancels when calculating the cor(i,j)
  for (int i=0; i<nnp; i++)
    for (int j=0; j<=i; j++)
      cov(i,j)= cov(j,i)= (sumij(i,j) - (sumi[i]*sumi[j])/sumw) / sumw;

  TH2D* hcor= new TH2D ("hcor", "Nuisance Parameter Correlation Matrix", nnp, 0.0, nnp, nnp, 0.0, nnp);
  hcor->SetMinimum(-1.0);
  hcor->SetMaximum( 1.0);
  hcor->GetXaxis()->SetLabelSize(0.013);
  hcor->GetYaxis()->SetLabelSize(0.013);
  hcor->GetZaxis()->SetLabelSize(0.025);

  int nv= nnp*(nnp-1)/2;
  vector<int> cori(nv), corj(nv);
  vector<Double_t> corv(nv);
  int k=0, np=0;
  for (int i=0; i<nnp; i++) {
    cor(i,i)= 1.0;
    hcor->SetBinContent(i+1,i+1,1.0);
    hcor->GetXaxis()->SetBinLabel(i+1,npnames[i].c_str());
    hcor->GetYaxis()->SetBinLabel(i+1,npnames[i].c_str());
    double fracerr= (cov(i,i)<=0.0) ? 0.0 : fabs((sumw*sqrt(cov(i,i)))/sumi[i]);
    if (verbose>=1)
      cout << npnames[i] << " mean = " << sumi[i]/sumw
           << " RMS = " << (cov(i,i)<=0.0 ? 0.0 : sqrt(cov(i,i))) << " fraction = " << fracerr << endl;
    if (fracerr<1e-6) cov(i,i)= 0.0;  // very small spread, could cause rounding errors
    for (int j=i+1; j<nnp; j++,k++) {
      Double_t Viijj= cov(i,i)*cov(j,j);
      if (Viijj>0.0) {
        Double_t rho= cov(i,j)/sqrt(Viijj);
        cor(i,j)= cor(j,i)= rho;
        corv[k]= fabs(rho);
        hcor->SetBinContent(i+1,j+1,rho);
        hcor->SetBinContent(j+1,i+1,rho);
      }
      cori[k]= i;
      corj[k]= j;
    }
    if (i<npoi) np=k+1;
  }
  assert(k==nv);

  if (hcor_ret) *hcor_ret= hcor;
  else delete hcor;

  // Pick most-correlated values
  vector<int> index(nv);
  TMath::Sort (np,    &corv.front(),    &index.front(),    true);
  TMath::Sort (nv-np, &corv.front()+np, &index.front()+np, true);

  if (maxCorrPOI>=0.0) {
    for (k=0; k<np; k++)
      if (corv[index[k]]<maxCorrPOI) break;
  } else
    k= np;
  cout << Form("npoi=%d, np=%d, nv=%d, k=%d",npoi,np,nv,k) << endl;
  for (int l=np; l<nv; l++, k++) {
    index[k]= np+index[l];
    if (maxCorr>=0.0 && corv[index[k]]<maxCorr) break;
  }
  nv= k;
  if (maxTop>=0 && nv>maxTop) nv= maxTop;

  topcor.resize(nv);
  itop.resize(nv);
  jtop.resize(nv);
  for (k=0; k<nv; k++) {
    int l=index[k], i=cori[l], j=corj[l];
    topcor[k]= cor(i,j);
    itop[k]= i;
    jtop[k]= j;
  }
  return nv;
}

void PValueCalculator::PlotPOIvsTS (TPad* canvas, Double_t tsdata, const char* poiname)
{
  TString ntcut= "!is_data";
  if (statsTree->GetBranch("is_alt")) ntcut += " && !is_alt";
  if (cut.Length())                   ntcut += Form(" && (%s)",     cut.Data());
  if (statsTree->GetBranch("weight")) ntcut  = Form("(%s)*weight",ntcut.Data());

  Long64_t nv= statsTree->Draw (poiname, ntcut, "groff");
  TH1* hy= statsTree->GetHistogram();
  if (nv==0 || !hy) return;
  Double_t ylo= hy->GetXaxis()->GetXmin(), yhi= hy->GetXaxis()->GetXmax();
  if (poimin[0]>ylo) ylo= poimin[0];
  if (poimax[0]<yhi) yhi= poimax[0];
  canvas->SetLogz();
  canvas->SetRightMargin(0.12);
  TH1D* hx= NewTSHist("hx","",tsdata);  // just for dimensions
  TH2D* hp= new TH2D ("hp", Form(";-2ln#Lambda;%s",poiname),
                      hx->GetNbinsX(), hx->GetXaxis()->GetXmin(), hx->GetXaxis()->GetXmax(),
                      100, ylo, yhi);
  delete hx;
  Project (statsTree, hp, Form("%s:2.0*(ml_cond-ml)",poiname), ntcut);
  IncludeOverflows(hp);
  hp->Draw("colz");
  TLine* line= 0;
  if (!isnan(tsdata) && tsdata<=hp->GetXaxis()->GetXmax()) {
    line= new TLine (tsdata, ylo, tsdata, yhi);
    line->SetLineWidth(3);
    line->Draw();
  }
  PrintPlot(canvas);
  canvas->SetRightMargin(0.03);
  canvas->SetLogz(0);
  delete line;
  delete hp;
}

void PValueCalculator::PlotTS2 (TPad* canvas, int tsType, Double_t tsdata, Double_t muhat,
                                const char* poiname, const char* hname, const char* hename)
{
  bool twoSided= (tsType==6);
  if (!isnan(tsdata)) {
    if      (tsType==6)                  tsdata= (muhat>=0 ? 1.0 : -1.0) * max(tsdata,0.0);
    else if (tsType==5 && muhat<bPOI[0]) tsdata= 0.0;
  }
  TString ntcut= "!is_data";
  if (statsTree->GetBranch("is_alt")) ntcut += " && !is_alt";
  if (cut.Length())                   ntcut += Form(" && (%s)",     cut.Data());
  if (statsTree->GetBranch("weight")) ntcut  = Form("(%s)*weight",ntcut.Data());

  const char* tsexp= tsType==6 ? Form("(%s>=0.0 ? 1.0 :-1.0) * 2.0*max(ml_cond-ml,0.0)",poiname) :
                     tsType==5 ? Form("(%s>=%g) ? 2.0*(ml_cond-ml) : -1e-14",poiname,bPOI[0])
                               :                 "2.0*(ml_cond-ml)";
  if (verbose>=0) cout << "Test statistic #"<<tsType<<" expression: "<<tsexp<<" cut: "<<ntcut<<endl;
  canvas->SetLogy();
  TH1D* h2= NewTSHist("h2", hname, tsdata, twoSided);
  Project (statsTree, h2, tsexp, ntcut);
  PlotTS (canvas, h2, 0, 0, tsType, tsdata);
  delete h2;

  canvas->SetLogy(0);
  TH1D* he2= NewTSHist("he2", hename, tsdata, twoSided, 1);
  Project (statsTree, he2, tsexp, ntcut);
  PlotTSdiff (canvas, he2, tsType, tsdata);
  delete he2;
}

void PValueCalculator::PlotVar (TPad* canvas, const char* name, const char* dname,
                                const char* varcut, Long64_t idata, Double_t genval)
{
  Long64_t nv= statsTree->Draw (name, varcut);
  TH1F* hn= dynamic_cast<TH1F*>(statsTree->GetHistogram());
  if (nv==0 || !hn) return;
  hn->SetTitle(name);
  Double_t xnlo=hn->GetXaxis()->GetXmin(), xnhi=hn->GetXaxis()->GetXmax();

  Double_t dval=NaN, derr=NaN;
  if (idata>=0) {
    dval= GetVal(statsTree,dname,idata);
    if (!isnan(dval))
      derr= GetVal(statsTree,Form("%s_err",dname),-1,-1.0);
  }

  vector<Double_t> statvals;
  vector<string> statnames;
  TLine* gline=0;
  if (!isnan(genval)) {
    statnames.push_back("Generated");
    statvals.push_back(genval);
    if (genval>=xnlo && genval<=xnhi) {
      gline= new TLine (genval, hn->GetMinimum(), genval, 0.8*hn->GetMaximum());
      gline->SetLineWidth(3);
      gline->Draw();
    }
  }
  if (!isnan(dval)) {
    statnames.push_back (!isnan(derr) && derr>0.0 ? "Data fit" : "Data");
    statvals.push_back(dval);
  }
  if (!isnan(derr) && derr>0.0) {
    statnames.push_back((optimize & kMinosData) ? "Minos err" : "Migrad err");
    statvals.push_back(derr);
  }
  if (!statnames.empty()) {
    TString func="0";
    for (int i=0, n=int(statnames.size()); i<n; i++)
      func += Form("+[%d]",i);
    TF1* statfun= new TF1 ("c", func, xnlo-1.0, xnlo-1.0);  // just for stats box
    for (size_t i=0, n=statnames.size(); i<n; i++) {
      statfun->SetParName (i, statnames[i].c_str());
      statfun->SetParameter (i, statvals[i]);
    }
    if (hn->GetListOfFunctions()) hn->GetListOfFunctions()->Add(statfun);
  }

  TLine* dline=0;
  if (!isnan(dval) && !isnan(derr) && derr>0.0) {
    Double_t xnb=(xnhi-xnlo)/hn->GetNbinsX(), sumw= hn->GetSum();
    if (verbose>=1) cout << name << ": sumw="<<sumw<<", mu="<<dval<<", err="<<derr<<", range="<<xnlo<<":"<<xnhi<<endl;
    TF1* gaus= new TF1 ("gaus", "gausn", xnlo, xnhi);
    gaus->SetParameters   (sumw*xnb,   dval,       derr);
    gaus->SetNpx(3*hn->GetNbinsX());
    gaus->SetLineColor(kRed);
    gaus->Draw("lsame");
  } else if (!isnan(dval) && dval>=xnlo && dval<=xnhi) {
    dline= new TLine (dval, hn->GetMinimum(), dval, 0.8*hn->GetMaximum());
    dline->SetLineWidth(3);
    if (gline) dline->SetLineStyle(2);
    dline->Draw();
  }

  PrintPlot (canvas);
  delete dline;
  delete gline;
}


void
PValueCalculator::PlotResult(TPad* canvas)
{
  TString massName;
  if (invMass>0.0) massName.Form(" at %g GeV",invMass);

  canvas->SetLogy();
  RooStats::SamplingDistribution* null= res->GetNullDistribution();
  if (!null) return;
  RooStats::SamplingDistribution* alt=  res->GetAltDistribution();

  bool twoSided= (testStatType==6);
  TString varName= null->GetVarName();
  TH1D* h= NewTSHist ("hnull", Form ("%s%s", varName.Data(), massName.Data()), 2.0*res->GetTestStatisticData(), twoSided);
  TH1D* ha= (alt && alt->GetSize()) ? dynamic_cast<TH1D*>(h->Clone("halt")) : 0;
  FillSamples (h, null, !twoSided);
  if (ha) FillSamples (ha, alt);
  TH1D* hb= 0;
  if (fitLeadbetter) {
    hb= NewTSHist ("hfit", Form ("%s%s", varName.Data(), massName.Data()), 2.0*res->GetTestStatisticData(), twoSided, 2);
    FillSamples (hb, null, !twoSided);
  }
  TF1* fit= 0;
  PlotTS (canvas, h, ha, hb, testStatType, 2.0*res->GetTestStatisticData(), &fit);
  delete hb;
  delete ha;
  if (!atlasStyle) canvas->SetLeftMargin(0.1);

  canvas->SetLogy(0);
  TH1D* he= NewTSHist ("htoys", Form ("Comparison with #chi^{2}%s", massName.Data()), 2.0*res->GetTestStatisticData(), twoSided, 1);
  FillSamples (he, null, !twoSided);  // rebinned for comparison
  PlotTSdiff (canvas, he, testStatType, 2.0*res->GetTestStatisticData(), fit);

  canvas->SetLogy();
  TH1D* hp= NewTSHist ("hpval", Form ("p-value%s", massName.Data()), 2.0*res->GetTestStatisticData(), twoSided);
  PlotPValue (canvas, hp, res, testStatType, fit);
  delete fit;
  canvas->SetLogy(0);
  canvas->SetRightMargin(0.07);

  if (statsTree && doStatsTree!=-2) {
    const Long64_t nent= statsTree->GetEntries();
    Int_t is_data=0;
    statsTree->SetBranchStatus("*",0);
    SetBranchAddress (statsTree, "is_data", &is_data);
    Long64_t idata;
    for (idata=0;idata<nent;idata++){
      statsTree->GetEntry(idata);
      if (is_data) break;
    }
    statsTree->SetBranchStatus("*",1);
    statsTree->GetEntry(idata);  // get rest of entry

    bool justCond= !statsTree->GetBranch("have_uncond") && statsTree->GetBranch("have_cond");
    bool haveBoth=  statsTree->GetBranch("have_uncond") && statsTree->GetBranch("have_cond");
    const string skip_var[] = {"weight", "is_data", "is_alt", "sample", "have_uncond", "have_cond"};
    const string* const skip_var_end = skip_var + sizeof(skip_var)/sizeof(skip_var[0]);
    const string   is_var[] = {"mu_init"};
    const string* const   is_var_end =   is_var + sizeof(  is_var)/sizeof(  is_var[0]);

    TIter ileaf(statsTree->GetListOfLeaves());
    vector<string> poinames, varnames, npnames;
    while (const TLeaf* leaf= dynamic_cast<const TLeaf*>(ileaf())) {
      const char* name= leaf->GetName();
      if (std::find (skip_var, skip_var_end, name) != skip_var_end) continue;
      TString sname= name;
      if (sname.EndsWith("_err"))                                   continue;
      if (sname.Contains("gamma_stat_"))                            continue;
      if (sname.BeginsWith("globObs_"))                             continue;
      bool have_err= statsTree->GetBranch(sname+"_err");
      bool is_cond= sname.EndsWith("_cond");
      if (have_err && (justCond ? !is_cond : is_cond))              continue;
      if (justCond) sname.Remove (sname.Length()-5);
      if      (!have_err || std::find (is_var, is_var_end, name) != is_var_end)
        varnames.push_back(sname.Data());
      else if (std::find (poiName.begin(), poiName.end(), sname.Data()) != poiName.end())
        poinames.push_back(sname.Data());
      else
        npnames.push_back(sname.Data());
    }
    if (poinames.empty()) {
      poinames.push_back("mu"); // probably won't work, but at least it'll give a better error message
      cerr << "Could not identify POI in ntuple - try "<<poinames[0]<<endl;
    }

    Double_t tsdata=NaN;
    vector<Double_t> muhat(poinames.size(),NaN), muerr(poinames.size(),NaN);
    if (is_data) {
      for (size_t i=0, n=poiName.size(); i<n; i++) {
        muhat[i]= GetVal (statsTree, poinames[i].c_str());
        muerr[i]= GetVal (statsTree, Form("%s_err",poinames[i].c_str()));
      }
      if (haveBoth) tsdata= GetVal (statsTree, "ml_cond") - GetVal (statsTree, "ml");
    } else
      idata= -1;
    statsTree->ResetBranchAddresses();
    bool haveWeight=statsTree->GetBranch("weight"), haveAlt=statsTree->GetBranch("is_alt");

    TString varcut= "!is_data";
    if (haveAlt)      varcut += " && !is_alt";
    if (cut.Length()) varcut += Form(" && (%s)",      cut.Data());
    if (haveWeight)   varcut  = Form("(%s)*weight",varcut.Data());

    if (haveBoth) {
      for (size_t i=0, n=poinames.size(); i<n; i++)
        PlotPOIvsTS (canvas, 2.0*tsdata, poinames[i].c_str());
      Long64_t nfail1= statsTree->GetEntries (Form("%s && ml_cond-ml<(%g)", varcut.Data(), negVal));
      Long64_t nfail2= 0;
      if (poinames.size()==1 && nullPOI.size()==1 && poimin[0]<nullPOI[0]) {
        nfail2= statsTree->GetEntries (Form("%s && %s>=0 && ml_cond-ml<(%g)", varcut.Data(), poinames[0].c_str(), negVal));
        if (testStatType!=5) PlotTS2 (canvas, 5, 2.0*tsdata, muhat[0], poinames[0].c_str(), h->GetTitle(), he->GetTitle());
        if (testStatType!=6) PlotTS2 (canvas, 6, 2.0*tsdata, muhat[0], poinames[0].c_str(), h->GetTitle(), he->GetTitle());
        if (testStatType!=2) PlotTS2 (canvas, 2, 2.0*tsdata, muhat[0], poinames[0].c_str(), h->GetTitle(), he->GetTitle());
      }
      if (verbose>=0) cout << nfail1 << " / " << nfail2 << " fit failures for all/+ve muhat toys" << endl;
    }

    canvas->SetLogy();
    gStyle->SetStatFontSize (.025);
    gStyle->SetStatFont (42);
    gStyle->SetStatX (0.99-canvas->GetRightMargin());
    gStyle->SetStatY (0.99-canvas->GetTopMargin());
    gStyle->SetStatW (0.10);
    gStyle->SetOptStat (1110);
    gStyle->SetOptFit (10001);

    for (size_t i=0, n=min(poiName.size(),bPOI.size()); i<n; i++) {
      PlotVar (canvas, poiName[i].c_str(),  poiName[i].c_str(),  varcut, idata, bPOI[i]);
    }
    for (size_t i=0; i<varnames.size(); i++) {
      PlotVar (canvas, varnames[i].c_str(), varnames[i].c_str(), varcut, idata);
    }

    gStyle->SetStatW (0.12);

    TString npcut= "!is_data";
    if (haveAlt)           npcut += " && !is_alt";
    if (plotMuHatPositive) {
      for (size_t i=0,n=poinames.size(); i<n; i++)
        npcut += Form(" && %s>=0.0",poinames[i].c_str());
    }
    if (cut.Length())      npcut += Form(" && (%s)",     cut.Data());
    if (haveWeight)        npcut  = Form("(%s)*weight",npcut.Data());
    if (verbose>=0) cout << "Ntuple selector: "<<npcut<<endl;
    for (size_t i=0; i<npnames.size(); i++) {
      PlotVar (canvas, npnames[i].c_str(), npnames[i].c_str(), npcut, idata);
    }
    canvas->SetLogy(0);
    gStyle->SetOptStat(0);

    // Find top pairs (|rho|>0.2) of correlated nuisance parameters (plus mu)
    npnames.insert(npnames.begin(), poinames.begin(), poinames.end());  // add POIs at the beginning
    vector<int> itop, jtop;
    vector<Double_t> topcor;
    TH2D* hcor=0;
    int nv= TopCorr (npnames, poinames.size(), itop, jtop, topcor, &hcor, 0.05, 0.2, 200);
    TAttPad save(*canvas);
    if (hcor) {
      canvas->SetBottomMargin(0.188);
      canvas->SetLeftMargin(0.184);
      canvas->SetRightMargin(0.1);
      hcor->Draw("colz");
      PrintPlot (canvas);
      delete hcor;
      save.Copy(*canvas);
    }
    if (nv>0) {
      gStyle->SetTitleYOffset (1.6);
      canvas->SetLeftMargin(0.13);
      canvas->SetRightMargin(0.12);
      canvas->SetLogz();
      for (int k=0;k<nv;k++) {
        int i=itop[k], j=jtop[k];
        if (verbose>=0)
          cout << "Nuisance parameters " << npnames[i] << ", " << npnames[j]
               << " have correlation coefficient " << topcor[k] << endl;
        Long64_t nr= statsTree->Draw (Form("%s:%s",npnames[i].c_str(),npnames[j].c_str()), npcut, "colz");
        TH1* hc= statsTree->GetHistogram();
        if (nr==0 || !hc) continue;
        hc->SetTitle(Form("%s vs %s (#rho=%.1f%%)", npnames[i].c_str(), npnames[j].c_str(), 100.0*topcor[k]));
        TLatex t;
        t.SetTextSize(0.03);
        t.SetNDC();
        t.DrawLatex(0.75,0.93,Form("#rho=%.1f%%",100.0*topcor[k]));
        PrintPlot (canvas);
      }
      save.Copy(*canvas);
      canvas->SetLogz(0);
    }
  }
  delete h;
  delete he;
}




//==============================================================================
// Constructors and destructor
//==============================================================================

PValueCalculator::PValueCalculator (const char* name)
  : WorkspaceCalculator (name)
{
  Init();
}

PValueCalculator::PValueCalculator (const char* name, int argc, const char* const* argv)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (argc, argv);
  cmdArgs= Join(vector<string>(argv+1,argv+argc)," ");
}

PValueCalculator::PValueCalculator (const char* name, const char* args)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (args);
  cmdArgs= args;
}

PValueCalculator::~PValueCalculator()
{
}

void PValueCalculator::Usage() const
{
  cout << "Usage:"
       << "\n  "<< GetName()<<" [ WORKSPACE-FILE.root | RESULT-FILES.root ] \\"
       << "\n      -w workspace|'combWS' -m modelConfig|'"<<modelSBName<<"' -d dataset|'"<<dataName<<"' -S Seed|FILE.root[:N] \\"
       << "\n      -p jobNum -t ntoys[:skipToys][,ntoysAlt]|"<<ntoys<<",0 -W nProofWorkers -M invMass -z poimin:poimax|"<<Join(poimin,poimax)<<" -Z vars|"<<Join(*varSettings,",",2)<<" \\"
       << "\n      -u muBkg[,muSig]|"<<Join(bPOI,":")<<","<<Join(sbPOI,":")<<" -O opt-level|";
  PrintOptimize(optimize,0);
  cout << " -A (asymptotic) -T testStatType|"<<testStatType<<" -e fitTol -H (plot) -c (ntuple cut) \\"
       << "\n      -n (no NP ntuple) -I #samples -v (verbose) -q (quiet) -r OUTPUT-RESULT-FILE"
       << endl;
}

int PValueCalculator::ParseArgs (int argc, const char* const* argv)
{
  vector<string> optimizeStrings;

  for (int i= 1; i<argc; i++) {
    if (argv[i][0]!='-') {
      fileName.push_back(argv[i]);
      continue;
    }
    for (const char* a= argv[i]+1; *a;) {
      switch (const char c= *a++) {
        case 'h': case '?':                    break;
        case 'H': plotResult++;                break;
        case 'A': calculatorType= (calculatorType!=2) ? 2 : 3;
                  break;
        case 'v': verbose++;                   break;
        case 'q': verbose--;                   break;
        case 'D': detailedOutput++;            break;
        case 'n': doStatsTree--;               break;
        case 'f': force = true;                break;
        default : if (!*a && i+1<argc) a= argv[++i];
        switch (c) {
        case 'w': wsName          = a;
                  altResultName   = true;      break;
        case 'm': modelSBName     =        a;  break;
        case 'd': dataName        =        a;  break;
        case 'r': resultFileName  =        a;  break;
        case 'S': seedName        =        a;  break;
        case 'c': cut             =        a;  break;
        case 'L': initialSnapshot =        a;
                  optimize |= kLoadInitialSnapshot;
                  break;
        case 'G': nullMLESnapshot =        a;
                  optimize |= kUseCondSnapshot;
                  break;
        case 'x': if (!ReadFile(a, wsEdit)) return 1;
                  break;
        case 'O': optimizeStrings.push_back(a);break;
        case 'C': a= Scan (a, poiName, ':');   break;
        case 'X': a= Scan (a, wsEdit, ';', true); break;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
        case 'l': initialSnapshotTS=       a;  break;
        case 'o': if (optimizeTS.Length()) optimizeTS += "+";
                  optimizeTS     +=        a;  break;
#endif
        default : const char* ai= a;
        switch (c) {
        case 'p': jobNum          = Strtol(a); break;
        case 't': if (*a==':') {
                    skipToys= Strtol(++a);
                    ntoys= skipToys+1;
                  } else {
                    ntoys= Strtol(a);
                    if (*a==':') skipToys= Strtol(++a);
                  }
                  if (*a==',') ntoysAlt= Strtol(++a);
                  nToysLimit      = true;      break;
        case 'u': a= Scan (a, bPOI, ':');
                  if (*a==',') a= Scan (a+1, sbPOI, ':');
                  break;
        case 'W': nworkers        = Strtol(a); break;
        case 'M': invMass         = Strtod(a); break;
        case '0': a= Scan (a, poimin);         break;  // obsolete
        case '1': a= Scan (a, poimax);         break;  // obsolete
        case 'z': a= Scan (a, poimin, poimax); break;
        case 'Z': a= Scan (a, *varSettings);   break;
        case 'T': testStatType    = Strtol(a); break;
        case 'e': if (*a=='p' && *++a) {
                    fitPrec= Strtod(a);
                    if (*a==',') a++;
                  }
                  if (*a && *a!=',') fitTol= Strtod(a); // Target EDM=0.001*fitTol. Default RooFit fitTol=1 (overrides Minuit fitTol=0.01).
                  if (*a==',' && *++a) {
                    fitTol0= fitTol;
                    fitTol= Strtod(a);
                  }
                  break;
        case 's': skipPlot        = Strtol(a); break;
        case '~': asimovSig       = Strtod(a); break;
        case 'I': {
          samplerType= 0;
          if (strcmp(a,"N")==0) {  // just null, but use IS code
            nStdDevOverlap= 0.0;
            a++;
          } else {
            double is             = Strtod(a);
            if (is>0.0 && is<5.0 && std::floor(is)!=is)
              nStdDevOverlap= is;
            else
              samplerType= TMath::Nint(is);
          }
          break;
        }
        case 'F': maxFunctionCalls= Strtol(a); break;
        case 'a': newParameterRanges= Strtod(a);
                  optimize |= kAdjustRanges;   break;
        case 'E': fixStatErrorMin = Strtod(a); 
                  optimize |= kFixStatError;   break;
        case 'i': if (*a==',')         interpCodeNorm= -1;
                  else                 interpCodeNorm=  Strtol(a);
                  if (*a==',' && *++a) interpCodeShape= Strtol(a);
                  else                 interpCodeShape= -1;
                  optimize |= kInterpolateErrors;
                  break;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
        case 'P': minosSetNum     = Strtol(a);
                  if (*a=='/') minosSetSize = Strtol(++a);
                  else         minosSetSize = 1;
                  if (*a) {
                    minosSetNum= minosSetSize= 0;
                    a= Scan (ai, minosSetNames);
                  }
                  optimize |= kMinosData;
                  break;
#endif
        default:          cerr << argv[0] << ": invalid option -"       << c              << endl; return 1; }
        if      (*a)    { cerr << argv[0] << ": invalid option value -" << c << " " << ai << endl; return 1; }
        else if (a==ai) { cerr << argv[0] << ": missing -" << c << " option value"        << endl; return 1; }}
        a= ""; // go to next arg
      }
    }
  }

  for (size_t i= 0, n=optimizeStrings.size(); i<n; i++)
    if (!SetOptimize(optimizeStrings[i].c_str(),"-O")) return 1;

  if (seedName.Length()>0) {
    char* a= 0;
    long int setSeed = strtol (seedName.Data(), (char**)&a, 10);
    if (*a == '\0') {
      seed= setSeed;
      seedName= "";
    }
  } else if (jobNum>=0)
    seed= -1;

  if (fileName.empty()) {
    SetDefaults();
    Usage();
    return 1;
  }

  return 0;
}
