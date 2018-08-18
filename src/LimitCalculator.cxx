/*==============================================================================
$Id: LimitCalculator.cxx 691299 2015-08-25 22:20:50Z adye $

Class to perform cross-section upper limit calculation using toys.
This calculator is used by StandardHypoTestInv.
Uses the WorkspaceCalculator base class for control and setup and
RooStats::HypoTestInverter (inverted hypothesis test) to perform
the scan and compute limits.

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

  useCLs           scan for CLs (otherwise for CLs+b)

  npoints:         number of points to scan , for autoscan set npoints = -1

  scanMin,scanMax: min/max value to scan in case of fixed scans
                   (if min >= max, try to find automatically)

  ntoys:           number of toys to use
  nToysRatio:      ratio of S+B/B toys (default is 2)

  plotResult:      plot result of tests at each point (TS distributions)
  writeResult:     write result of scan (default is true)

Author: Tim Adye, based on $ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C
                  from ROOT trunk revision 41812 (5.32.00-rc1+).

==============================================================================*/

#include "Utilities/WorkspaceCalculatorConfig.h"
#include "LimitCalculator.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TList.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStopwatch.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"

#include "RooStats/HypoTestCalculatorGeneric.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::isnan;

const int                nlim = 5;
const Double_t    lim   [nlim]= { -2.0,       -1.0,       0.0,       1.0,        2.0     };
const char* const limnam[nlim]= {"-2 sigma", "-1 sigma", "median", "+1 sigma", "+2 sigma"};


void LimitCalculator::Init()
{
  // Don't use conditional snapshot since each point needs different profiled parameters.
  // Don't ever skip the data fit because HypoTestInverter::RunOnePoint rejects results with p-value=NaN.
  optimize          &= ~(kUseCondSnapshot|kSkipDataTS);
  jobSet             = 0;
  bPOIdefault        = 0.0;
  useCLs             = true;
  npoints            = 50;    // Number of mu points to scan in [scanMin:scanMax]
  scanMin            = 0;
  scanMax            = 10;
  confidenceLevel    = 0.95;
  rebuild            = false; // re-do extra toys for computing expected limits and rebuild test stat distributions
                              // N.B this requires much more CPU (factor is equivalent to nToyToRebuild)
  nToyToRebuild      = 100;   // number of toys used to rebuild
  nToysRatio         = 2.0;   // ratio Ntoys S+B/B
  firstPoint         = 0;
  lastPoint          = -1;    // Set >=0 to scan only a subset of the points
  expCLsFromSB       = true;  // expected CLs limits written to ntuple from CLs+b, else from background
  resultFileName     = "Results.root";
  resultName         = "result_mu";
  showAtlasLabel     = false;
  runLimit           = true;
  runToys            = true;

  res                = 0;
  InitCalculator();
}

void LimitCalculator::InitCalculator()
{
  limitCalc = 0;
}

void LimitCalculator::DeleteCalculator()
{
  delete limitCalc;
  WorkspaceCalculator::DeleteCalculator();
}

void LimitCalculator::SetDefaults()
{
  if (testStatType<0) testStatType= 3;  // test statistic type (3=one-sided)
  if (res) res->UseCLs(useCLs);
  WorkspaceCalculator::SetDefaults();
}


void LimitCalculator::ShowParms() const
{
     cout << GetName()
          << " wsName=\""         << wsName      << "\""
          << ", modelSBName=\""   << modelSBName << "\""
          << ", modelBName=\""    << modelBName  << "\""
          << ", dataName=\""      << dataName    << "\"";
     PrintOptimize (optimize & ~(kAdjustRanges|kFixStatError), ", ");
     if (optimize & kAdjustRanges)   cout << ", AdjustRanges "<<newParameterRanges<<" sigma";
     if ((optimize & kFixStatError) && fixStatErrorMin>=0)
                                     cout << ", FixStatError "<<100.0*fixStatErrorMin<<"%";
     if (detailedOutput)             cout << ", DetailedOutput";
     if (dropNegative)               cout << ", dropNegative";
     cout << ", npoints="         << npoints;
     if (!pointsToScan.empty()) cout << ", pointsToScan=" << Join(pointsToScan);
     if      (lastPoint==firstPoint)
     cout << ", scanPoint="       << firstPoint;
     else if (lastPoint> firstPoint)
     cout << ", scanPoints="      << firstPoint << '-' << lastPoint;
     cout << ", ntoys="           << ntoys;
     if (skipToys>0)
     cout << ", skipToys="        << skipToys;
     if (ntoysAlt>0)
     cout << ", ntoysAlt="        << ntoysAlt;
     cout << ", " << (useCLs ? "CLs" : "CLs+b");
     if (nworkers>0)
     cout << ", nworkers="        << nworkers;
     if (invMass>0.0)
     cout << ", invMass="         << invMass;
     if (!poimin.empty())
     cout << ", poimin="          << Join(poimin);
     if (!poimax.empty())
     cout << ", poimax="          << Join(poimax);
     if (scanMax>scanMin)
     cout << ", scanMin="         << scanMin
          << ", scanMax="         << scanMax;
     cout << ", confidenceLevel=" << confidenceLevel;
     cout << endl;
}


int LimitCalculator::RunOverWorkspace()
{
  // run the calculator
  int     ok= SetupWorkspace();
  if (ok) ok= SetupMinimizer();
  if (ok) ok= SetupModel();
  if (ok) ok= SetupInitialParms();
  if (ok) ok= SetupTestStat();
  if (ok) ok= SetupSampler();
  if (ok) ok= SetupCalculator();
  if (ok) ok= SetupScan();
  if (ok) ok= RunCalculatorProof();
  if (ok) ok= RebuildExpected();
  if (ok) ok= GetTestStatInfo();
  DeleteCalculator();   // don't use from here on
  return ok;
}


int LimitCalculator::SetupScan()
{
   limitCalc = new RooStats::HypoTestInverter(*calc);
   limitCalc->SetConfidenceLevel(confidenceLevel);
   limitCalc->UseCLs(useCLs);
   limitCalc->SetVerbose(verbose+1);

   if (npoints > 0) {
      if (scanMin >= scanMax) {
         // if no min/max given scan between MLE and +4 sigma
         scanMin = poihat;
         scanMax = poihat +  4 * poierr;
      }
      if (lastPoint<firstPoint) lastPoint = npoints-1;
      if (lastPoint==firstPoint) {
        Double_t poival;
        if (!pointsToScan.empty()) poival = pointsToScan[firstPoint];
        else                       poival = scanMin + (firstPoint+1)*(scanMax-scanMin)/npoints;
        cout << "Doing a single point at "<<poival<<endl;
        limitCalc->SetFixedScan(1,poival,poival);
      } else if (!pointsToScan.empty()) {
        cout << "Doing a fixed "<<lastPoint-firstPoint+1<<"-point scan over points ";
        for (Int_t i=firstPoint; i<=lastPoint; i++) {
          if (i>firstPoint) cout << ',';
          cout << pointsToScan[i];
        }
        cout << endl;
      } else {
        Int_t npt = lastPoint-firstPoint+1;
        Double_t step= (scanMax-scanMin)/npoints;
        Double_t lo= scanMin + (firstPoint+1)*step;
        Double_t hi= scanMin + ( lastPoint+1)*step;
        cout << "Doing a fixed "<<npt<<"-point scan in interval : " << lo << " , " << hi << endl;
        limitCalc->SetFixedScan(npt,lo,hi);
      }
   }
   else {
      //poi->setMax(10*int( (poihat+ 10 *poierr )/10 ) );
     RooRealVar* poi = dynamic_cast<RooRealVar*>(nullModel->GetParametersOfInterest()->first());
     cout << "Doing an  automatic scan  in interval : " << poi->getMin() << " , " << poi->getMax() << endl;
   }
   return 1;
}

int LimitCalculator::RunCalculator()
{
  TStopwatch tw;

  if (lastPoint>firstPoint && !pointsToScan.empty()) {
    for (Int_t i=firstPoint; i<=lastPoint; i++) {
      if (!limitCalc->RunOnePoint (pointsToScan[i])) {
        cout << "Loop interupted because of failed status" << endl;
        return 0;
      }
    }
    cout << "Now retrieve result"<<endl;
  }

  res = limitCalc->GetInterval();
  res->SetName(resultName);   // use a fixed name to simplify reading back. Name is "result_mu" for compatibility.
  if (verbose>=0) {
    cout << "Time to perform limit scan: "; tw.Print();
  }
  return (res ? 1 : 0);
}

int LimitCalculator::RebuildExpected()
{
  if (!rebuild) return 1;
  TStopwatch tw;
  RooStats::SamplingDistribution * limDist = limitCalc->GetUpperLimitDistribution (true, nToyToRebuild);
  if (verbose>=0) {
    cout << "Time to rebuild distributions: "; tw.Print();
  }
  if (limDist) {
    cout << "Rebuilt expected limits: ";
    for (int i=0; i<nlim; i++) {
      if (i>0) cout << ", ";
      cout << limnam[i] << ": " << limDist->InverseCDF(ROOT::Math::normal_cdf(lim[i]));
    }
    cout << endl;
  } else
    cout << "ERROR : failed to re-build distributions" << endl;
  delete res;
  res = limitCalc->GetInterval();
  return (res ? 1 : 0);
}


void LimitCalculator::AdjustResults (bool resultWasRead)
{
  WorkspaceCalculator::AdjustResults (resultWasRead);
  if (res) saved.Add (res);
}


//==============================================================================
// Other methods
//==============================================================================

TTree* LimitCalculator::GetLimits()
{
  if (!res) return 0;

  Int_t nEntries = res->ArraySize();
  if (nEntries<=0) return 0;

  bool asymptotic= !res->GetNullTestStatDist(0) && !res->GetAltTestStatDist(0);

  // sort the values in x
  vector<Double_t> xvals(nEntries);
  for (int i= 0; i<nEntries; i++) xvals[i]= res->GetXValue(i);
  vector<int> index(nEntries);
  TMath::Sort (nEntries, &xvals.front(), &index.front(), false);

  Double_t resMin=0, resMax=0;
  Int_t nToysSig=0, nToysBkg=0;
  Int_t nobs=0, nexp=0;
  if (asymptotic) {
    if (verbose>=0) cout << "Asymptotic points:";
    for (int j= 0; j<nEntries; j++) {
      int i= index[j];
      Double_t x= xvals[i];
      if (j==0 || x<resMin) resMin= x;
      if (j==0 || x>resMax) resMax= x;
      if (res->GetResult(i)->CLb()>0.0) nobs++;
      if (verbose>=0) cout << " " << x;
    }
    if (verbose>=0) cout << endl;
    nexp= nobs;
  } else {
    nToysSig=nEntries/2; nToysBkg=nEntries/2;  // round average
    if (verbose>=0) cout << "Points:toysSB/toysB: ";
    for (int j= 0; j<nEntries; j++) {
      int i= index[j];
      Double_t x= xvals[i];
      if (j==0 || x<resMin) resMin= x;
      if (j==0 || x>resMax) resMax= x;
      Int_t nToysNull= res->GetNullTestStatDist(i) ? res->GetNullTestStatDist(i)->GetSize() : 0;
      Int_t nToysAlt=  res->GetAltTestStatDist (i) ? res->GetAltTestStatDist (i)->GetSize() : 0;
      nToysSig += nToysNull;
      nToysBkg += nToysAlt;
      if (nToysNull>0 &&  nToysAlt>0) nexp++;
      if (nToysNull>0 && (nToysAlt>0 || !useCLs) && res->GetResult(j)->HasTestStatisticData()) nobs++;
      if (verbose>=0 && j>0) {
        if (j%5==0) cout << endl << "                     ";
        else        cout << ' ';
      }
      if (verbose>=0) cout << x << ':' << nToysNull << '/' << nToysAlt;
    }
    if (verbose>=0) cout << endl;
    nToysSig /= nEntries;  // average
    nToysBkg /= nEntries;  // average
  }
    
  if (verbose>=0) {
    RooArgSet* args= res->GetParameters();
    const RooRealVar* poi = dynamic_cast<RooRealVar*>(args->first());
    cout << nEntries << ' ' << poi->GetName()
         << " points in range " << resMin << " to " << resMax
         << " (fit range " << poi->getMin() << " to " << poi->getMax() << "), "
         << nToysSig << " null toys, "
         << nToysBkg << " alt toys per point" << endl;
    delete args;
  }

  Double_t obsUL=0.0, obsUL_err=0.0, exp[nlim];
  Double_t obserr=0.0, expserr[nlim], expberr[nlim];
  for (int i=0; i<nlim; i++) {
    exp[i]= expserr[i]= expberr[i]= 0.0;
  }

  if (!haveWeights && !asymptotic)
    ResampleErrors (res, (nobs>=2 ? &obserr : 0), (nexp>=2&&useCLs ? expserr : 0), (nexp>=2 ? expberr : 0), nToysResample);

  if (nobs>=2) {
    obsUL=        res->UpperLimit();
    obsUL_err=    res->UpperLimitEstimatedError();
    cout << Form("Observed upper    limit = %-9g +/- %-9g +/- %-9g", obsUL, obsUL_err, obserr) << endl;
  }

  if (nexp>=2) {
    for (int i=0; i<nlim; i++) {
      Double_t expb= res->GetExpectedUpperLimit(lim[i]);
      if (useCLs) {
        res->UseCLs(false);
        double cl= 1.0 - (1.0-confidenceLevel)*RooStats::SignificanceToPValue(lim[i]);
        res->SetConfidenceLevel (cl);
        exp[i]= res->GetExpectedUpperLimit(0);
        res->SetConfidenceLevel(confidenceLevel);
        res->UseCLs(true);
        cout << Form("expected %-8s limit = %-9g +/- %-9g from CLs+b at %.1f%%, limit from bkg = %-9g +/- %g",
                     limnam[i], exp[i], expserr[i], 100.0*cl, expb, expberr[i]) << endl;
      } else {
        cout << Form("expected %-8s limit = %-9g +/- %g", limnam[i], expb, expberr[i]) << endl;
      }
      if (!expCLsFromSB || !useCLs) {
        exp[i]= expb;
        expserr[i]= expberr[i];
      }
    }
    if (useCLs && verbose>=0) {
      if (expCLsFromSB) cout << "expected limits from CLs+b written to ntuple" << endl;
      else              cout << "expected limits from bkg written to ntuple" << endl;
    }
  }

  // Rongkun added those parts to directly write out the result in a histogram
  if( !writeResult )
  {
      TFile* fileOut = TFile::Open (resultFileName.c_str(), "recreate");
      fileOut->cd();
      TH1D h_lim("limit","limit",7,0,7);
      h_lim.SetBinContent(1,res->UpperLimit());
      h_lim.SetBinContent(2,res->GetExpectedUpperLimit(0));
      h_lim.SetBinContent(3,res->GetExpectedUpperLimit(+2));
      h_lim.SetBinContent(4,res->GetExpectedUpperLimit(+1));
      h_lim.SetBinContent(5,res->GetExpectedUpperLimit(-1));
      h_lim.SetBinContent(6,res->GetExpectedUpperLimit(-2));

      h_lim.GetXaxis()->SetBinLabel(1, "Observed");
      h_lim.GetXaxis()->SetBinLabel(2, "Expected");
      h_lim.GetXaxis()->SetBinLabel(3, "+2sigma");
      h_lim.GetXaxis()->SetBinLabel(4, "+1sigma");
      h_lim.GetXaxis()->SetBinLabel(5, "-1sigma");
      h_lim.GetXaxis()->SetBinLabel(6, "-2sigma");
      h_lim.GetXaxis()->SetBinLabel(7, "Global status");

      fileOut->Write();
      fileOut->Close();
      //delete fileOut;
  }
  //hist.Fill(6,res->GetExpectedUpperLimit());
  //Global Satus

  UInt_t initialSeed= seed;
  TTree* bandTree= new TTree ("band", "StandardHypoTestInv");
  bandTree->SetDirectory(0);   // memory resident ntuple
  bandTree->Branch ("invMass",      &invMass);
  bandTree->Branch ("band2sigDown", &exp[0]);
  bandTree->Branch ("band1sigDown", &exp[1]);
  bandTree->Branch ("bandMedian",   &exp[2]);
  bandTree->Branch ("band1sigUp",   &exp[3]);
  bandTree->Branch ("band2sigUp",   &exp[4]);
  bandTree->Branch ("obsUL",        &obsUL);
  bandTree->Branch ("obsUL_err",    &obserr);
  bandTree->Branch ("band2sigDown_err", &expserr[0]);
  bandTree->Branch ("band1sigDown_err", &expserr[1]);
  bandTree->Branch ("bandMedian_err",   &expserr[2]);
  bandTree->Branch ("band1sigUp_err",   &expserr[3]);
  bandTree->Branch ("band2sigUp_err",   &expserr[4]);
  bandTree->Branch ("muhat",        &poihat);
  bandTree->Branch ("muhat_err",    &poierr);
  bandTree->Branch ("muMin",        &resMin);
  bandTree->Branch ("muMax",        &resMax);
  bandTree->Branch ("muBkg",        &bPOI);
  bandTree->Branch ("confidenceLevel",&confidenceLevel);
  bandTree->Branch ("nPointsToScan",&nEntries);
  bandTree->Branch ("nToysSig",     &nToysSig);
  bandTree->Branch ("nToysBkg",     &nToysBkg);
  bandTree->Branch ("useCLs",       &useCLs);
  bandTree->Branch ("wsfile",       &wsfile);
  bandTree->Branch ("optimize",     &optimize);
  bandTree->Branch ("seed",         &initialSeed);
  bandTree->Branch ("args",         &cmdArgs);
  bandTree->Branch ("poiName",      &poiName);
  bandTree->Branch ("tstype",       &testStatType);
  bandTree->Fill();   // just one entry so we can combine later
  bandTree->ResetBranchAddresses();

  saved.Add (bandTree);
  return bandTree;
}


int LimitCalculator::ReadResult (TFile* f, const char* objname)
{
  RooStats::HypoTestInverterResult* radd = 0;
  f->GetObject (objname, radd);
  if (!radd) return 1;
  if (verbose>=1) PrintResult (radd);
  if (!res) {
    RooArgSet* p= radd->GetParameters();
    RooRealVar* v= dynamic_cast<RooRealVar*>(p->first());
    delete p;
    if (!v) {
      cerr << "Bad HypoTestInverterResult::"<<objname<<" object from file " << f->GetName() << ": no scanned parameter" << endl;
      delete radd;
      return 2;
    }
    res= new RooStats::HypoTestInverterResult (radd->GetName(), *v, radd->ConfidenceLevel());
    res->SetTitle(radd->GetTitle());
    res->UseCLs(useCLs);
  }
  if (verbose<1) RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  bool ok= res->Add(*radd);
  if (verbose<1) RooMsgService::instance().getStream(1).addTopic(RooFit::Eval);
  if (!ok) {
    cerr << "Could not merge HypoTestInverterResult::"<<objname<<" object from file " << f->GetName() << endl;
    delete radd;
    return 2;
  }
  delete radd;
  return 0;
}


int LimitCalculator::ReadTree (TFile* f, UInt_t& thisSeed)
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

  Double_t invMass1=0.0, poihat1=0.0, poierr1=0.0;
  vector<Double_t> bPOI1;
  vector<string>* poiName1=0;
  Int_t testStatType1=-1;
  string* wsfile1=0;
  ULong64_t optimize1=0;
  bandTree->SetBranchAddress ("invMass", &invMass1);
  if (bandTree->GetBranch("muhat"))     bandTree->SetBranchAddress ("muhat",     &poihat1);
  if (bandTree->GetBranch("muhat_err")) bandTree->SetBranchAddress ("muhat_err", &poierr1);
  if (bandTree->GetBranch("tstype"))    bandTree->SetBranchAddress ("tstype",    &testStatType1);
  if (bandTree->GetBranch("wsfile"))    bandTree->SetBranchAddress ("wsfile",    &wsfile1);
  SetBranchAddress (bandTree, "muBkg", &bPOI1, true, NaN);
  if (bandTree->GetBranch("optimize"))  bandTree->SetBranchAddress ("optimize",  &optimize1);
  if (bandTree->GetBranch("seed"))      bandTree->SetBranchAddress ("seed",      &thisSeed);
  if (bandTree->GetBranch("poiName"))   bandTree->SetBranchAddress ("poiName",   &poiName1);
  bool haveVec= (bPOI1.size()==1);
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
      if (!force) return 4;
    }
  }
  if (! bPOI1.empty() && !isnan( bPOI1[0])) {
    if ( bPOI.empty())  bPOI=  bPOI1;
    if ( bPOI1!= bPOI && !(haveVec &&  bPOI.size()>1 &&  bPOI1.size()==1 &&  bPOI1[0]== bPOI[0])) {
      cerr << "results with different Bkg POI settings: " << Join( bPOI,":")<< " and " << Join( bPOI1,":") << endl;
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


Double_t LimitCalculator::ExpectedPValue (RooStats::HypoTestResult* r, Double_t* experr, Double_t* q0exp, double nsig)
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


void
LimitCalculator::ResampleErrors (const RooStats::HypoTestInverterResult* rold,
                                 Double_t* obserr, Double_t* expserr, Double_t* expberr,
                                 int nToys)
{
  if (!(obserr || expserr || expberr)) return;
  Double_t obs=0.0, obs2=0.0, expb[nlim], expb2[nlim], exps[nlim], exps2[nlim];
  for (int i=0; i<nlim; i++) {
    exps[i]= exps2[i]= expb[i]= expb2[i]= 0.0;
    if (obserr) *obserr=     0.0;
    if (expserr) expserr[i]= 0.0;
    if (expberr) expberr[i]= 0.0;
  }
  if (nToys<2) return;
  cout << "StandardHypoTestInv: resample ";
  if (obserr)  cout << "observed, ";
  if (expserr) cout << "expected S+B, ";
  if (expberr) cout << "expected bkg, ";
  cout << "sampling distributions " << nToys << " times to estimate errors" << endl;
  RooFit::MsgLevel errlevel= RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  for (int j=0; j<nToys; j++) {
    RooStats::HypoTestInverterResult* r= ResampleToy(rold);
    if (obserr) {
      Double_t v= r->UpperLimit();
      obs  += v;
      obs2 += v*v;
    }
    if (expberr) {
      for (int i=0; i<nlim; i++) {
        Double_t v= r->GetExpectedUpperLimit(lim[i]);
        expb [i] += v;
        expb2[i] += v*v;
      }
    }
    if (expserr) {   // only for CLs
      r->UseCLs(false);
      for (int i=0; i<nlim; i++) {
        double cl= 1.0 - (1.0-confidenceLevel)*RooStats::SignificanceToPValue(lim[i]);
        r->SetConfidenceLevel (cl);
        Double_t v= r->GetExpectedUpperLimit(0);
        exps [i] += v;
        exps2[i] += v*v;
      }
      r->UseCLs(true);
      r->SetConfidenceLevel(confidenceLevel);
    }
    delete r;
  }
  RooMsgService::instance().setGlobalKillBelow(errlevel);
  if   (obserr)     *obserr= sqrt (fabs (obs2     - (obs    *obs    )/nToys) / (nToys-1));
  for (int i=0; i<nlim; i++) {
    if (expserr) expserr[i]= sqrt (fabs (exps2[i] - (exps[i]*exps[i])/nToys) / (nToys-1));
    if (expberr) expberr[i]= sqrt (fabs (expb2[i] - (expb[i]*expb[i])/nToys) / (nToys-1));
  }
}

RooStats::HypoTestInverterResult*
LimitCalculator::ResampleToy (const RooStats::HypoTestInverterResult* rold)
{
  RooStats::HypoTestInverterResult* r= dynamic_cast<RooStats::HypoTestInverterResult*>(rold->Clone());
  const int nEntries = r->ArraySize();
  for (int i= 0; i<nEntries; i++) {
    RooStats::SamplingDistribution* null= ResampleToy (r->GetNullTestStatDist(i));
    if (null) r->GetResult(i)->SetNullDistribution(null);
    RooStats::SamplingDistribution* alt=  ResampleToy (r->GetAltTestStatDist(i));
    if (alt)  r->GetResult(i)->SetAltDistribution(alt);
  }
  return r;
}


int LimitCalculator::LimitSamples()
{
  int npt= 0;
  const int nEntries = res->ArraySize();
  for (int i= 0; i<nEntries; i++) {
    RooStats::SamplingDistribution* null= LimitSamples (res->GetNullTestStatDist(i), ntoys,    jobNum);
    if (null) res->GetResult(i)->SetNullDistribution(null);
    RooStats::SamplingDistribution* alt=  LimitSamples (res->GetAltTestStatDist(i),  ntoysAlt, jobNum);
    if (alt)  res->GetResult(i)->SetAltDistribution(alt);
    if (null || alt) npt++;
  }
  if (npt>0) cout << "Reduced to "<<ntoys<<" S+B toys and "<<ntoysAlt<<" bkg toys for "<<npt<<" points"<<endl;
  return npt;
}


void
LimitCalculator::PrintResult (RooStats::HypoTestInverterResult* r)
{
  RooArgSet* args= r->GetParameters();
  const RooRealVar* poi = dynamic_cast<RooRealVar*>(args->first());

  cout << "============== ";
  r->Print();
  const int nEntries = r->ArraySize();
  cout << "Result";
  if (!wsfile.empty()) cout << " for workspace file " << wsfile;
  if (invMass>0.0) cout << " at " << invMass << " GeV";
  if (!poiName.empty()) cout << " " << Join(poiName);
  else                  cout << " mu";
  if (poierr!=0.0) cout << " fitted = " << poihat << " +/- " << poierr;
  cout << " has "<<nEntries<<" points"<<endl;

  // sort the values in x
  vector<Double_t> xvals(nEntries);
  for (int j= 0; j<nEntries; j++) xvals[j]= r->GetXValue(j);
  vector<int> index(nEntries);
  TMath::Sort (nEntries, &xvals.front(), &index.front(), false);
  for (int j= 0; j<nEntries; j++) {
    int i= index[j];
    cout << "Point "<<j+1<<"/"<<nEntries<<" at "<<poi->GetName()<<'='<<r->GetXValue(i)<<":";
    r->GetResult(i)->Print();

    RooStats::SamplingDistribution* s = r->GetExpectedPValueDist(i);
    if (s) {
      const vector<double>& values = s->GetSamplingDistribution();
      if (values.size()>=2) {
        double p[nlim], q[nlim];
        for (int k=0; k<nlim; k++) {
          p[k]= ROOT::Math::normal_cdf(lim[k]);
          q[k]= 0;
        }
        double* x = const_cast<double*>(&values[0]);
        TMath::Quantiles (values.size(), nlim, x,q,p,false);
        cout << "Expected limits: ";
        for (int k=0; k<nlim; k++) {
          if (k>0) cout << ", ";
          cout << limnam[k] << ": " << q[k];
        }
        cout << endl;
      }
      PrintSamplingDistribution (s, "p-value");
      delete s;
    }
    PrintSamplingDistribution (r->GetNullTestStatDist(i), "null   ");
    PrintSamplingDistribution (r->GetAltTestStatDist(i),  "alt    ");
    cout << endl;
  }
  delete args;
}


int LimitCalculator::DropBadSamples()
{
  haveWeights= 0;
  int npt= 0;
  const int nEntries = res->ArraySize();
  for (int i= 0; i<nEntries; i++) {
    RooStats::SamplingDistribution* null= DropBadSamples (res->GetNullTestStatDist(i));
    if (null) res->GetResult(i)->SetNullDistribution(null);
    RooStats::SamplingDistribution* alt=  DropBadSamples (res->GetAltTestStatDist(i));
    if (alt)  res->GetResult(i)->SetAltDistribution(alt);
    if (null || alt) npt++;
  }
  if (npt>0) cout << "Dropped bad samples from "<<npt<<" points"<<endl;
  return npt;
}

int LimitCalculator::ApplyCuts()
{
  haveWeights= 0;
  int npt= 0;
  const int nEntries = res->ArraySize();
  for (int i= 0; i<nEntries; i++) {
    RooStats::SamplingDistribution* null= ApplyCuts (res->GetNullTestStatDist(i));
    if (null) res->GetResult(i)->SetNullDistribution(null);
    RooStats::SamplingDistribution* alt=  ApplyCuts (res->GetAltTestStatDist(i));
    if (alt)  res->GetResult(i)->SetAltDistribution(alt);
    if (null || alt) npt++;
  }
  if (npt>0) cout << "Applied cuts to "<<npt<<" points"<<endl;
  return npt;
}


void LimitCalculator::PlotResult(TPad* canvas)
{
  gStyle->SetOptTitle(1);
  canvas->SetTopMargin(0.07);
  canvas->SetLeftMargin(0.09);

  TString poiname= "mu";
  if (RooArgSet* args= res->GetParameters()) {
    if (const RooRealVar* poi = dynamic_cast<RooRealVar*>(args->first()))
      poiname= poi->GetName();
    delete args;
  }

  const int nEntries = res->ArraySize();
  if (nEntries<=0) return;

  TString massName;
  if (invMass>0.0) massName.Form(" at %g GeV",invMass);

  bool pvok= true;
  for (int i= 0; i<nEntries; i++) {
    if (!isnan(res->GetResult(i)->NullPValue()) &&
        !isnan(res->GetResult(i)->AlternatePValue()) &&
        res->GetResult(i)->CLb()>0.0)
      continue;
    pvok= false;
    break;
  }

  if (pvok) {
    RooStats::HypoTestInverterPlot *plot = new RooStats::HypoTestInverterPlot(res);
    TGraphErrors* graph= plot->MakePlot();
    // force p-value axis range 0..1 (prevent CLs=-1 from screwing it up)
    TH1* h= graph->GetHistogram();
    h->SetMinimum(0.0);
    h->SetMaximum(1.05);
    h->SetTitle(Form("%s%s;%s;p value",plot->GetTitle(),massName.Data(),poiname.Data()));
    h->Draw();
    plot->Draw("2CL CLb SAME");
    Double_t sz= 1.0 - res->ConfidenceLevel();
    TLine line (h->GetXaxis()->GetXmin(), sz, h->GetXaxis()->GetXmax(), sz);
    line.SetLineColor(kRed);
    line.Draw();
    PrintPlot (canvas);
    delete graph;

    if (useCLs) {
      canvas->SetLogy();
      res->UseCLs(false);
      TMultiGraph* gexp = plot->MakeExpectedPlot();
      gexp->Draw("A");
      gexp->GetXaxis()->SetTitle(poiname);
      gexp->GetYaxis()->SetTitle("Expected CL_{s+b}");

      double x1 = gexp->GetXaxis()->GetXmin();
      double x2 = gexp->GetXaxis()->GetXmax();
      for (int i=0; i<nlim; i++) {
        double lsz= sz * RooStats::SignificanceToPValue(lim[i]);
        TLine* lline = new TLine (x1, lsz, x2, lsz);
        lline->SetLineColor(kRed);
        lline->Draw();
      }
      PrintPlot (canvas);
      canvas->SetLogy(0);
      delete gexp;
      res->UseCLs(true);
    }
    delete plot;
  }

  // sort the values in x
  vector<Double_t> xvals(nEntries);
  for (int j= 0; j<nEntries; j++) xvals[j]= res->GetXValue(j);
  vector<int> index(nEntries);
  TMath::Sort (nEntries, &xvals.front(), &index.front(), false);

  Double_t maxts= 0;
  TVectorD poix(nEntries), tsy(nEntries), nullp(nEntries), altp(nEntries);
  int nvec= 0;
  for (int j=0; j<nEntries; j++) {
    int i= index[j];
    const RooStats::HypoTestResult* r= res->GetResult(i);
    if (!r->HasTestStatisticData()) continue;
    Double_t ts= 2.0*r->GetTestStatisticData();
    poix[nvec]= res->GetXValue(i);
    tsy[nvec]= ts;
    nullp[nvec]= r->NullPValue();
    altp[nvec]=  r->AlternatePValue();
    if (j==0 || ts>maxts) maxts= ts;
    nvec++;
  }
  if (nvec>0) {
    TString tsvar=tsName();
    poix.ResizeTo(nvec);
    tsy.ResizeTo(nvec);
    nullp.ResizeTo(nvec);
    altp.ResizeTo(nvec);
    int sides= (testStatType==2) ? 2 : 1;
    {
      TGraph g (poix, tsy);
      g.SetTitle(Form("%s%s;%s;%s",tsvar.Data(),massName.Data(),poiname.Data(),tsvar.Data()));
      g.SetLineColor(1);
      g.SetLineWidth(2);
      g.SetMarkerStyle(20);
      g.SetMarkerColor(1);
      g.Draw("alp");
      PrintPlot (canvas);
    }
    {
      TGraph g (poix, nullp);
      g.SetTitle(Form("Null p-value%s;%s;p-value",massName.Data(),poiname.Data()));
      g.SetLineColor(1);
      g.SetLineWidth(2);
      g.SetMarkerStyle(20);
      g.SetMarkerColor(1);
      g.Draw("alp");
      PrintPlot (canvas);
      TGraph gi (nullp, poix);  // TGraph to do interpolation
      cout << "Null "<<sides<<"-sided limits: ";
      for (int i=0; i<=3; i++) {
        if (i>0) cout << ", ";
        cout << i << "sigma" << " " << gi.Eval(sides*ROOT::Math::normal_cdf(-i));
      }
      cout << endl;
    }
    {
      TGraph g (poix, altp);
      g.SetTitle(Form("Alt p-value%s;%s;p-value",massName.Data(),poiname.Data()));
      g.SetLineColor(1);
      g.SetLineWidth(2);
      g.SetMarkerStyle(20);
      g.SetMarkerColor(1);
      g.Draw("alp");
      PrintPlot (canvas);
      TGraph gi (altp, poix);  // TGraph to do interpolation
      cout << "Alt  "<<sides<<"-sided limits: ";
      for (int i=0; i<=3; i++) {
        if (i>0) cout << ", ";
        cout << i << "sigma" << " " << gi.Eval(sides*ROOT::Math::normal_cdf(-i));
      }
      cout << endl;
    }
  }

  canvas->SetLogy();

  Double_t xhi= 1.5*maxts;
  if (xhi<5) xhi= 5;
  for (int j=0; j<nEntries; j++) {
    int i= index[j];
    const RooStats::HypoTestResult* r= res->GetResult(i);
    RooStats::SamplingDistribution* null= r->GetNullDistribution();
    RooStats::SamplingDistribution* alt=  r->GetAltDistribution();
    if (!null && !alt) continue;

    bool twoSided= (testStatType==6);
    TString varName= null->GetVarName();
    TH1D* h= NewTSHist ("hnull", Form ("%s for %s=%g%s", varName.Data(), poiname.Data(), res->GetXValue(i), massName.Data()),
                        2.0*r->GetTestStatisticData(), twoSided);
    TH1D* ha= (alt && alt->GetSize()) ? dynamic_cast<TH1D*>(h->Clone("halt")) : 0;
            FillSamples (h, null, !twoSided);
    if (ha) FillSamples (ha, alt, !twoSided);
    PlotTS (canvas, h, ha, 0, testStatType, 2.0*r->GetTestStatisticData());
  }

  for (int j=0; j<nEntries; j++) {
    int i= index[j];
    RooStats::SamplingDistribution* pvals= res->GetExpectedPValueDist(i);
    TH1F hp ("pvalue", Form("expected p-value for %s=%g%s;p-value", poiname.Data(), res->GetXValue(i), massName.Data()), 100, 0, 1);
    FillSamples (&hp, pvals, kBlack);
    hp.Draw();
    PrintPlot (canvas);
    delete pvals;
  }
}




//==============================================================================
// Constructors and destructor
//==============================================================================

LimitCalculator::LimitCalculator (const char* name)
  : WorkspaceCalculator (name)
{
  Init();
}

LimitCalculator::LimitCalculator (const char* name, int argc, const char* const* argv)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (argc, argv);
  cout << "HAHAHA "<<  writeResult << endl;
  cmdArgs= Join(vector<string>(argv+1,argv+argc)," ");
}

LimitCalculator::LimitCalculator (const char* name, const char* args)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (args);
  cmdArgs= args;
}

LimitCalculator::~LimitCalculator()
{
}

void LimitCalculator::Usage() const
{
  cout << "Usage:"
       << "\n  "<< GetName()<<" [ WORKSPACE-FILE.root | RESULT-FILES.root ] \\"
       << "\n      -w workspace|'combWS' -m modelConfig|'"<<modelSBName<<"' -d dataset|'"<<dataName<<"' -S Seed|FILE.root[:N] \\"
       << "\n      -a (all points) -p jobNum -n nPointsToScan|"<<npoints<<" -N points/job -P pointsToScan -t ntoys[:skipToys][,ntoysAlt]|"<<ntoys<<",0 \\"
       << "\n      -W nProofWorkers -M invMass -0 scanMin|"<<scanMin<<" -1 scanMax|"<<scanMax<<" -z poimin:poimax|"<<Join(poimin,poimax)<<" -Z vars|"<<Join(*varSettings,",",2)<<" -c confidenceLevel|"<<confidenceLevel<<" \\"
       << "\n      -B (use CLs+b, not CLs) -u muBkg|"<<Join(bPOI)<<" -O opt-level|";
  PrintOptimize(optimize,0);
  cout << " -A (asymptotic) -T testStatType|"<<testStatType<<" -e fitTol -H (plot) \\"
       << "\n      -v (verbose) -q (quiet) -r OUTPUT-RESULT-FILE"
       << endl;
}

int LimitCalculator::ParseArgs (int argc, const char* const* argv)
{
  const char*  pointsToScanArg = 0;
  double       scanLim         = 1e30;
  bool         all             = false;
  int          nPointsNow      = -1;
  vector<string> optimizeStrings;

  for (int i= 1; i<argc; i++) {
    if (argv[i][0]!='-') {
      fileName.push_back(argv[i]);
      continue;
    }
    for (const char* a= argv[i]+1; *a;) {
      switch (const char c= *a++) {
        case 'h': case '?':                    break;
        case 'a': all             = true;      break;
        case 'B': useCLs          = false;     break;
        case 'H': plotResult++;                break;
        case 'A': calculatorType  = 2;         break;
        case 'v': verbose++;                   break;
        case 'q': verbose--;                   break;
        case 'D': detailedOutput++;            break;
        case 'R': doStatsTree--;               break;
        case 'f': force = true;                break;
        // Rongkun add this line to save the writing memory..
        case 'K': writeResult        = false; break;
        default : if (!*a && i+1<argc) a= argv[++i];
        switch (c) {
        //those are for when there's a parameter following the -flag
        case 'w': wsName          = a;
                  altResultName   = true;      break;
        case 'm': modelSBName     =        a;  break;
        case 'd': dataName        =        a;  break;
        case 'r': resultFileName  =        a;  break;
        case 'P': pointsToScanArg =        a;  break;
        case 'S': seedName        =        a;  break;
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
        case 'n': npoints         = Strtol(a); break;
        case 'N': nPointsNow      = Strtol(a);
                  all             = false;     break;
        case 't': if (*a==':') {
                    skipToys= Strtol(++a);
                    ntoys= skipToys+1;
                  } else {
                    ntoys= Strtol(a);
                    if (*a==':') skipToys= Strtol(++a);
                  }
                  if (*a==',') ntoysAlt= Strtol(++a);
                  nToysLimit      = true;      break;
        case 'u': a= Scan (a, bPOI, ':');      break;
        case 'W': nworkers        = Strtol(a); break;
        case 'M': invMass         = Strtod(a); break;
        case '0': scanMin         = Strtod(a); break;
        case '1': scanMax=scanLim = Strtod(a); break;
        case '2': a= Scan (a, poimax);         break;  // obsolete
        case 'z': a= Scan (a, poimin, poimax); break;
        case 'Z': a= Scan (a, *varSettings);   break;
        case 'c': confidenceLevel = Strtod(a); break;
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
        case 'F': maxFunctionCalls= Strtol(a); break;
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

  if (pointsToScanArg) {  // -P 1,2,3,4,10/2,/10 -> 1,2,3,4,6,8,10,20,30,...,scanMax
    Double_t lastp= scanMin;
    for (char* next=0;; pointsToScanArg= next+1) {
      Double_t pointVal= strtod (pointsToScanArg, &next);
      if (*next=='/') {
        if (next==pointsToScanArg) pointVal= scanMax;
        pointsToScanArg= next+1;
        Double_t step= strtod (pointsToScanArg, &next);
        if (next>pointsToScanArg && step > 0.0) {
          Double_t hi= pointVal + 1e-6*step;
          Double_t p= lastp + step;
          for (; p < hi; p += step)
            if (p>=scanMin && p<scanLim+1e-6*step) pointsToScan.push_back(p<scanLim?p:scanLim);
        }
      } else if (pointVal>=scanMin && pointVal<=scanLim)
        pointsToScan.push_back(pointVal);
      if (*next != ',') break;
      lastp= pointVal;
    }
    npoints= pointsToScan.size();
    scanMax= pointsToScan[npoints-1];
  }

  if (jobNum<0 && nPointsNow<=0) all = true;
  if (jobNum<0) jobNum = 0;
  jobSet= jobNum;
  if (!all) {
    if (nPointsNow<=1) {
      firstPoint = lastPoint = jobNum % npoints;
      jobSet                 = jobNum / npoints;
    } else {
      jobNum *= nPointsNow;
      int np = ((npoints+nPointsNow-1)/nPointsNow)*nPointsNow;
      firstPoint = jobNum % np;
      jobSet     = jobNum / np;
      lastPoint = firstPoint+nPointsNow-1;
      if (lastPoint>=npoints) lastPoint = npoints-1;
    }
  }

  if (fileName.empty()) {
    SetDefaults();
    Usage();
    return 1;
  }

  return 0;
}
