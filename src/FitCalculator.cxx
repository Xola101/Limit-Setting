/*==============================================================================
$Id: FitCalculator.cxx 731458 2016-03-21 18:46:55Z adye $

Class to perform a simple fit to the data in a workspace.
Uses the WorkspaceCalculator base class for control and setup.

Main parameters:

  fileName         workspace file or list of result files
  wsName           workspace name
  modelSBName      ModelConfig object name
  dataName         observed dataset

  testStatType   = 0 Simple Likelihood Ratio (LEP)
                 = 1 Ratio of Profile Likelihood (Tevatron)
                 = 2 Profile Likelihood Ratio two-sided
                 = 3 Profile Likelihood Ratio one-sided           (i.e.   0  if mu_hat >= mu)
                 = 4 Maximum Likelihood Estimate (mu_hat)
                 = 5 Profile Likelihood Ratio one-sided discovery (i.e.   0  if mu_hat <= mu)
                 = 6 Profile Likelihood Ratio signed discovery    (i.e. -q0  if mu_hat <  mu)
                 = 7 Profile Likelihood Ratio signed              (i.e. -qmu if mu_hat >  mu)

Author: Tim Adye.

==============================================================================*/

#include "Utilities/WorkspaceCalculatorConfig.h"
#include "FitCalculator.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TList.h"
#include "TF1.h"
#include "TPad.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooFitResult.h"

#include "RooStats/RooStatsUtils.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/TestStatistic.h"

#ifdef USE_ProfileLikelihoodTestStatEnhanced
#include "ProfileLikelihoodTestStatEnhanced.h"
#endif


using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

void FitCalculator::Init() {
  optimize           &= ~kSkipDataTS;
  seed               = -1;    // don't use random numbers, so no point in randomising and saving seeds
  resultFileName     = "fit.root";
  resultName         = "impact_*,contrib_*";
  doStatsTree        = 0;
  nll                = NaN;
  errorAnalysis      = 0;
  nParmsNow          = 1;
  optionalResult     = true;
  nsigma             = 1.0;
  postfit            = true;
}


void FitCalculator::SetDefaults() {
  if (testStatType<0) testStatType= 2; // test statistic type (2=two-sided)
  if (errorAnalysis) optimizeTS= "F"+optimizeTS; // keep fit result
  if (verbose>=0) optimizeTS= "P"+optimizeTS;    // print fit result
  WorkspaceCalculator::SetDefaults();
}


void FitCalculator::ShowParms() const
{
     cout << GetName()
          << " wsName=\""         << wsName      << "\""
          << ", modelSBName=\""   << modelSBName << "\""
          << ", modelBName=\""    << modelBName  << "\""
          << ", dataName=\""      << dataName    << "\"";
     PrintOptimize (optimize & ~(kAdjustRanges|kMinosData), ", ");
     if (optimize & kAdjustRanges)   cout << ", AdjustRanges "<<newParameterRanges<<" sigma";
     if (optimize & kMinosData) {
                                     cout << ", MinosData";
       if (minosSetSize>0)           cout << " (set " << minosSetNum << "/" << minosSetSize << ")";
     }
     if (detailedOutput)             cout << ", DetailedOutput";
     if (invMass>0.0)
     cout << ", invMass="         << invMass;
     if (!poimin.empty())
     cout << ", poimin="          << Join(poimin);
     if (!poimax.empty())
     cout << ", poimax="          << Join(poimax);
     cout << endl;
}


int FitCalculator::RunOverWorkspace()
{
  // run the calculator
  int     ok= SetupWorkspace();
  if (ok) ok= SetupMinimizer();
  if (ok) ok= SetupModel();
  if (ok) ok= SetupInitialParms();
  if (ok) ok= SetupTestStat();
  if (ok) ok= RunCalculatorProof();
  if (ok) ok= GetTestStatInfo();
  if (!plotResult) DeleteCalculator();   // don't use from here on
  return ok;
}


int FitCalculator::RunCalculator()
{
   // run fit
   TStopwatch tw;

   bool dofit0= true;
   const RooArgSet* pois= nullModel->GetParametersOfInterest();
   const RooArgSet* impVars= minosSet ? minosSet : minosParams;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
   RooStats::ProfileLikelihoodTestStatEnhanced* ts= dynamic_cast<RooStats::ProfileLikelihoodTestStatEnhanced*>(testStat);
   if ((optimize & kMinosData) && ts) {
     if        (errorAnalysis==1) {
       if (!postfit) {
         RooArgSet* poiMinos= dynamic_cast<RooArgSet*>(impVars->selectCommon(*pois));
         if (poiMinos->getSize() > 0) {
           ts->UseMinos(poiMinos);
           kept.Add(poiMinos);
         } else {
           delete poiMinos;
           dofit0= false;
         }
       }
     } else if (errorAnalysis==2)
       ts->UseMinos(pois);
   }
#endif

   RooArgSet* snap0= 0;
   if (errorAnalysis) snap0= dynamic_cast<RooArgSet*>(nullParams->snapshot());

   RooArgSet poiSet= *nullSnapshot;

   if (optimize & kAltMinos) {

     minos_results= impVars->snapshot();
     minos_results->setName("minos_results");
     RooAbsCollection* fitted_poi= pois->snapshot();
     RooArgList fitparms;   // use list so can hold different snapshots of the same var if POI included in minosSet
     fitparms.add(*fitted_poi);
     fitparms.add(*minos_results);

     nll= DoFit (poiSet, 1, &fitparms);

     cout << "    ----> Run findSigma for ";
     minos_results->Print();

     TList impactList;
     RooAbsCollection* fitted_lo= pois->snapshot();
     RooAbsCollection* fitted_hi= pois->snapshot();
     for (RooFIter it= minos_results->fwdIterator(); RooAbsArg* a= it.next();) {
       RooRealVar* p= dynamic_cast<RooRealVar*>(a);
       if (!p) continue;
       *fitted_lo = *pois;
       *fitted_hi = *pois;

       findSigma (nll, p, nsigma, fitted_lo, fitted_hi, "-PU", -1.0);

       if (errorAnalysis != 1) continue;

       if (verbose>=0) {
         cout << Form("At %s = %.3f (MLE%+.3f): ", a->GetName(), p->getVal()+p->getAsymErrorLo(), p->getAsymErrorLo());
         fitted_lo->Print("v");
         cout << Form("At %s = %.3f (MLE%+.3f): ", a->GetName(), p->getVal()+p->getAsymErrorHi(), p->getAsymErrorHi());
         fitted_hi->Print("v");
       }

       for (RooFIter pit= fitted_poi->fwdIterator(); RooAbsArg* vp= pit.next();) {
         RooRealVar* poi= dynamic_cast<RooRealVar*>(vp);
         if (!poi) continue;
         TString impactName= Form("impact_%s",poi->GetName());
         RooArgSet* impact= dynamic_cast<RooArgSet*>(impactList.FindObject(impactName));
         if (!impact) impactList.Add ((impact= new RooArgSet (impactName)));
         if (strcmp (poi->GetName(), p->GetName()) == 0) {
           impact->addOwned (*RooArgSet(*p).snapshot());
           continue;
         }
         Double_t poiVal= poi->getVal();
         Double_t poiLo=0.0, poiHi=0.0;
         if (const RooRealVar* v= dynamic_cast<const RooRealVar*>(fitted_lo->find(*poi)))
           poiLo= v->getVal() - poiVal;
         if (const RooRealVar* v= dynamic_cast<const RooRealVar*>(fitted_hi->find(*poi)))
           poiHi= v->getVal() - poiVal;
         cout << Form("%s +/-1sigma: %s = %.3f %+.3f %+.3f",a->GetName(),poi->GetName(),poiVal,poiHi,poiLo) << endl;

         RooRealVar* impVar= new RooRealVar (a->GetName(), Form("%s impact on %s/%s %s", a->GetName(), (poiHi>=0?"+":"-"), (poiLo>=0?"+":"-"), poi->GetName()),
                                             poiVal, poi->getMin(), poi->getMax());
         if (poiHi>=0.0)
           impVar->setAsymError (-fabs(poiLo), fabs(poiHi));
         else {
           impVar->setAsymError (-fabs(poiHi), fabs(poiLo));
           impVar->setAttribute ("anticorrelated");
         }
         if ((poiLo?poiLo:1.0)*(poiHi?poiHi:1.0) > 0.0)
           impVar->setAttribute ("same-sign");
         impact->addOwned (*impVar);
       }
     }
     delete fitted_lo;
     delete fitted_hi;
     if (verbose>=0)
       for (TIter iit= &impactList; const TObject* impact= iit();) {
         cout << "Impact on " << impact->GetName()+7 << ":" << endl;
         impact->Print("v");
       }
     saved.AddAll(&impactList);
     saved.Add(minos_results);

   } else {

     if (dofit0) nll= testStat->Evaluate (*data, poiSet);

#ifdef USE_ProfileLikelihoodTestStatEnhanced
   if (errorAnalysis==1 && ts) {
     bool resetNP= (ts->SetOptimize() & RooStats::ProfileLikelihoodTestStatEnhanced::kResetNP);

     ts->UseMinos(0,false);
     ts->SetOptimize("C-HUP",0);
//   ts->SetPrintLevel(verbose-1);

     RooArgSet* snap= dynamic_cast<RooArgSet*>(nullParams->snapshot());

     Double_t nup= nsigma*nsigma;

     for (RooFIter pit= pois->fwdIterator(); RooAbsArg* vp= pit.next();) {
       RooRealVar* poi= dynamic_cast<RooRealVar*>(vp);
       if (!poi) continue;
       Double_t poiVal= poi->getVal();
       RooArgSet* impact= new RooArgSet (Form("impact_%s",poi->GetName()));

       for (RooFIter it= impVars->fwdIterator(); const RooAbsArg* a= it.next();) {
         if (strcmp (poi->GetName(), a->GetName()) == 0) {
           impact->addOwned (*RooArgSet(*a).snapshot());
           continue;
         }
         const RooRealVar* p= dynamic_cast<const RooRealVar*>(a);
         if (!p) continue;
         RooArgSet* psnap= dynamic_cast<RooArgSet*>(RooArgSet(*p).snapshot());
         RooRealVar* pvar= dynamic_cast<RooRealVar*>(psnap->first());

         Double_t pval= 0.0, poi0= poiVal, nll0= nll;
         if (postfit)
           pval= p->getVal();
         else {
           if (resetNP) *nullParams= *snap0;
           std::cout  << " XXX I found you!" << std::endl;
           pvar->setVal (pval);
           nll0= ts->Evaluate (*data, *psnap);
           poi0= poi->getVal();
           cout << Form("%s = %g: %s = %g", p->GetName(), pval, poi->GetName(), poi0) << endl;
           *nullParams= *snap;
         }

         Double_t errLo= (postfit ? nsigma*p->getErrorLo() : -nsigma), valLo= pval + errLo;
         if (valLo <= p->getMin()) {
           cout << "Warning: MINOS -ve error is at limit: ";
           p->Print();
         }
         if (resetNP) *nullParams= *snap0;
         pvar->setVal (valLo);
         Double_t nllLo= ts->Evaluate (*data, *psnap);
         Double_t poiLo= poi->getVal() - poi0;
         Double_t dNllLo= 2.0*(nllLo-nll0);
         cout << Form("%s %+g = %g: %s = %g, delta = %+g, -2lnL = %.4f",
                      p->GetName(), errLo, valLo, poi->GetName(), poi->getVal(), poiLo, dNllLo);
         if (fabs(dNllLo/nup-1.0) > 0.05) cout << " - not "<<nup<<"UP!";
         cout << endl;
         *nullParams= *snap;

         Double_t errHi= (postfit ? nsigma*p->getErrorHi() :  nsigma), valHi= pval + errHi;
         if (valHi >= p->getMax()) {
           cout << "Warning: MINOS +ve error is at limit: ";
           p->Print();
         }
         if (resetNP) *nullParams= *snap0;
         pvar->setVal (valHi);
         Double_t nllHi= ts->Evaluate (*data, *psnap);
         Double_t poiHi= poi->getVal() - poi0;
         Double_t dNllHi= 2.0*(nllHi-nll0);
         if (verbose>=0)
           cout << Form("%s %+g = %g: %s = %g, delta = %+g, -2lnL = %.4f",
                        p->GetName(), errHi, valHi, poi->GetName(), poi->getVal(), poiHi, dNllHi);
         if (fabs(dNllHi/nup-1.0) > 0.05) cout << " - not "<<nup<<"UP!";
         cout << endl;
         *nullParams= *snap;

         RooRealVar* impVar= new RooRealVar (p->GetName(), Form("%s impact on %s/%s %s", p->GetName(), (poiHi>=0?"+":"-"), (poiLo>=0?"+":"-"), poi->GetName()),
                                             poi0, poi->getMin(), poi->getMax());
         if (poiHi>=0.0)
           impVar->setAsymError (-fabs(poiLo), fabs(poiHi));
         else {
           impVar->setAsymError (-fabs(poiHi), fabs(poiLo));
           impVar->setAttribute ("anticorrelated");
         }
         if ((poiLo?poiLo:1.0)*(poiHi?poiHi:1.0) > 0.0)
           impVar->setAttribute ("same-sign");
         impact->addOwned (*impVar);
         delete psnap;
       }
       if (verbose>=0) {
         cout << "Impact on " << poi->GetName() << ":" << endl;
         impact->Print("v");
       }
       saved.Add(impact);
     }

     delete snap;


   } else if (errorAnalysis==2 && ts) {
     bool resetNP= (ts->SetOptimize() & RooStats::ProfileLikelihoodTestStatEnhanced::kResetNP);

     ts->SetOptimize("-P",0);
//   ts->SetPrintLevel(verbose-1);

     RooArgSet* snap= dynamic_cast<RooArgSet*>(nullParams->snapshot());
     const RooArgList poiList= *pois;
     const Int_t npois= poiList.getSize();

     if (resetNP) *nullParams= *snap0;
     nll= testStat->Evaluate (*data, poiSet);

     RooArgList* poiSnap= dynamic_cast<RooArgList*>(poiList.snapshot());
     vector<RooArgSet*> contrib(npois);
     for (Int_t i= 0; i<npois; i++) {
       contrib[i]= new RooArgSet (Form("contrib_%s",poiList[i].GetName()));
       contrib[i]->addOwned (*poiSnap->at(i));
       saved.Add(contrib[i]);
     }
     poiSnap->releaseOwnership();
     delete poiSnap;

     for (RooFIter it= impVars->fwdIterator(); RooAbsArg* a= it.next();) {
       if (npois==1 && strcmp (poiList[0].GetName(), a->GetName()) == 0) continue;
       RooRealVar* p= dynamic_cast<RooRealVar*>(a);
       if (!p) continue;
       p->setConstant();
       if (verbose>=1) { cout << "    ----> Fix "; p->Print(); }
       nll= testStat->Evaluate (*data, poiSet);

       poiSnap= dynamic_cast<RooArgList*>(poiList.snapshot());
       for (Int_t i= 0; i<npois; i++) {
         RooRealVar* poi= dynamic_cast<RooRealVar*>(poiSnap->at(i));
         if (!poi) continue;
         if (strcmp (poi->GetName(), p->GetName()) == 0) continue;
         if (verbose>=0) {
           RooRealVar* pp= dynamic_cast<RooRealVar*>(snap->find(poiList[i]));
           if (!pp) pp= poi;
           cout << Form ("Fix %s: %s -> %g %+g %+g, change %+g %+g %+g",
                         p->GetName(), poi->GetName(),
                         poi->getVal(), poi->getAsymErrorLo(), poi->getAsymErrorHi(),
                         poi->getVal()-pp->getVal(),
                         pp->getAsymErrorLo()-poi->getAsymErrorLo(),
                         poi->getAsymErrorHi()-pp->getAsymErrorHi()) << endl;
         }
         RooAbsArg* poiCopy= dynamic_cast<RooAbsArg*>(poi->clone(p->GetName()));
         poiCopy->SetTitle(Form("%s removed from %s",p->GetName(),poi->GetName()));
         contrib[i]->addOwned (*poiCopy);
       }
       delete poiSnap;
       *nullParams= *snap;

     }

     delete snap;

     if (verbose>=0)
       for (Int_t i= 0; i<npois; i++) {
         cout << "POI " << poiList[i].GetName() << " with each parameter fixed:" << endl;
         contrib[i]->Print("v");
       }
   }
   }
#endif

   delete snap0;

   if (verbose>=0) {
     cout << "Time for fit: "; tw.Print();
   }
   return 1;
}



//==============================================================================
// Other methods
//==============================================================================

TTree* FitCalculator::GetLimits()
{
  int npoi= bPOI.size();
  bool twoSided= (testStatType==2 && !(npoi==1 && poimin[0]>=bPOI[0])); // really two-sided
  double sides= twoSided ? 2.0 : 1.0;

  double q= 2.0*nll, p, Z;
  if (npoi==1) {  // shortcut
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

  cout << "Asymptotic p-value " << p << ", significance " << Z << " ("<<npoi<<"DoF, "<<sides<<"-sided)" << endl;

  if (optimize & kSaveWorkspace) return 0;

  UInt_t initialSeed= seed;
  TTree* bandTree= new TTree ("band", GetName());
  bandTree->SetDirectory(0);   // memory resident ntuple
  bandTree->Branch ("tsdata",       &nll);
  bandTree->Branch ("invMass",      &invMass);
  bandTree->Branch ("muhat",        &poihat);
  bandTree->Branch ("muhat_err",    &poierr);
  bandTree->Branch ("muMin",        &poimin);
  bandTree->Branch ("muMax",        &poimax);
  bandTree->Branch ("tstype",       &testStatType);
  bandTree->Branch ("wsfile",       &wsfile);
  bandTree->Branch ("optimize",     &optimize);
  bandTree->Branch ("seed",         &initialSeed);
  bandTree->Branch ("args",         &cmdArgs);
  bandTree->Fill();   // just one entry so we can combine later
  bandTree->ResetBranchAddresses();

  saved.Add (bandTree);
  return bandTree;
}


int FitCalculator::ReadResult (TFile* f, const char* objname)
{
  RooAbsCollection* res = 0;
  f->GetObject (objname, res);
  if (!res) return 1;
  if (res->getSize() == 0) return 0;
  if (verbose>=1) cout << "File " << f->GetName() << " contains " << res->ClassName() << "::" << res->GetName() << " with " << res->getSize() << " entries" << endl;
  RooAbsCollection* r= dynamic_cast<RooAbsCollection*>(saved.FindObject(objname));
  if (!r) {
    r= new RooArgList (objname);
    saved.Add(r);
  }
  r->addClone(*res);
  delete res;
  return 0;
}


int FitCalculator::ReadTree (TFile* f, UInt_t& thisSeed)
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
  Int_t testStatType1=-1;
  string* wsfile1=0;
  ULong64_t optimize1=0;
  bandTree->SetBranchAddress ("invMass", &invMass1);
  if (bandTree->GetBranch("muhat"))     bandTree->SetBranchAddress ("muhat",     &poihat1);
  if (bandTree->GetBranch("muhat_err")) bandTree->SetBranchAddress ("muhat_err", &poierr1);
  if (bandTree->GetBranch("tstype"))    bandTree->SetBranchAddress ("tstype",    &testStatType1);
  if (bandTree->GetBranch("wsfile"))    bandTree->SetBranchAddress ("wsfile",    &wsfile1);
  if (bandTree->GetBranch("optimize"))  bandTree->SetBranchAddress ("optimize",  &optimize1);
  bandTree->GetEntry(0);
  if (wsfile1 && !wsfile1->empty()) {
    if      (wsfile.empty())
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
  if (poierr==0.0 && poierr1!=0.0) {
    poihat= poihat1;
    poierr= poierr1;
  }
  optimize |= optimize1;

  bandTree->ResetBranchAddresses();
  delete bandTree;
  return 0;
}


void
FitCalculator::PrintResult()
{
  cout << "Result";
  if (!wsfile.empty()) cout << " for workspace file " << wsfile;
  if (invMass>0.0) cout << " at " << invMass << " GeV";
  if (poierr!=0.0) cout << " with muhat = " << poihat << " +/- " << poierr;
  if (testStatType>=0) cout << " and test stat #" << testStatType;
  cout <<endl;
}


void FitCalculator::PlotResult (TPad* canvas)
{
  if (!nullPdf || !data) return;
  canvas->SetLeftMargin(0.1);
  canvas->SetRightMargin(0.05);
  PlotDataSet (canvas, *nullPdf, *data);
  DeleteCalculator();
}


//==============================================================================
// Constructors and destructor
//==============================================================================

FitCalculator::FitCalculator (const char* name)
  : WorkspaceCalculator (name)
{
  Init();
}

FitCalculator::FitCalculator (const char* name, int argc, const char* const* argv)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (argc, argv);
  cmdArgs= Join(vector<string>(argv+1,argv+argc)," ");
}

FitCalculator::FitCalculator (const char* name, const char* args)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (args);
  cmdArgs= args;
}

FitCalculator::~FitCalculator()
{
}


void FitCalculator::Usage() const
{
  cout << "Usage:"
       << "\n  "<< GetName()<<" [ WORKSPACE-FILE.root | RESULT-FILES.root ] \\"
       << "\n      -w workspace|'combWS' -m modelConfig|'"<<modelSBName<<"' -d dataset|'"<<dataName<<"' \\"
       << "\n      -M invMass -z poimin:poimax|"<<Join(poimin,poimax)<<" -C poiName|"<<Join(poiName,":")<<" \\"
       << "\n      -O opt-level|";
  PrintOptimize(optimize,0);
  cout
       <<" -T testStatType|"<<testStatType<<" -e fitTol \\"
       << "\n      -v (verbose) -q (quiet) -r OUTPUT-RESULT-FILE"
       << endl;
}

int FitCalculator::ParseArgs (int argc, const char* const* argv)
{
  bool doCond= false;
  vector<string> optimizeStrings;

  for (int i= 1; i<argc; i++) {
    if (argv[i][0]!='-') {
      fileName.push_back(argv[i]);
      continue;
    }
    for (const char* a= argv[i]+1; *a;) {
      switch (const char c= *a++) {
        case 'h': case '?':                    break;
        case 'v': verbose++;                   break;
        case 'q': verbose--;                   break;
        case 'H': plotResult++;                break;
        case 'D': detailedOutput++;            break;
        case 't': doStatsTree++;               break;
        case 'f': force = true;                break;
        case 'W': optimize |= kSaveWorkspace;  break;
        case 'A': errorAnalysis++;             break;
        case 'c': constrainedOnly = true;      break;
        default : if (!*a && i+1<argc) a= argv[++i];
        switch (c) {
        case 'w': wsName          = a;
                  altResultName   = true;      break;
        case 'm': modelSBName     =        a;  break;
        case 'd': dataName        =        a;  break;
        case 'r': resultFileName  =        a;  break;
        case 'L': initialSnapshot =        a;
                  optimize |= kLoadInitialSnapshot;
                  break;
        case 'x': if (!ReadFile(a, wsEdit)) return 1;
                  break;
        case 'O': optimizeStrings.push_back(a);break;
        case 'C': a= Scan (a, poiName, ':');   break;
        case 'X': a= Scan (a, wsEdit, ';', true); break;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
        case 'P': a= Scan (a, minosSetNames);  break;
        case 'l': initialSnapshotTS=       a;  break;
        case 'o': if (optimizeTS.Length()) optimizeTS += "+";
                  optimizeTS     +=        a;  break;
#endif
        default : const char* ai= a;
        switch (c) {
        case 'p': jobNum          = Strtol(a); break;
        case 'N': nParmsNow       = Strtol(a); break;
        case 'M': invMass         = Strtod(a); break;
        case 'u': a= Scan (a, bPOI, ':');
                  if (*a==',') a= Scan (a+1, sbPOI, ':');
                  doCond= true;
                  break;
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
        case 's': nsigma= Strtod(a);
                  if ((postfit= (*a=='s'))) a++;
                  break;
        default:          cerr << argv[0] << ": invalid option -"       << c              << endl; return 1; }
        if      (*a)    { cerr << argv[0] << ": invalid option value -" << c << " " << ai << endl; return 1; }
        else if (a==ai) { cerr << argv[0] << ": missing -" << c << " option value"        << endl; return 1; }}
        a= ""; // go to next arg
      }
    }
  }

  if (!postfit) constrainedOnly= true;

  if (jobNum >= 0) {
    minosSetNum=  jobNum;
    minosSetSize= nParmsNow;
  }

  for (size_t i= 0, n=optimizeStrings.size(); i<n; i++)
    if (!SetOptimize(optimizeStrings[i].c_str(),"-O")) return 1;

  if ((errorAnalysis || !minosSetNames.empty() || minosSetSize>0) && !(optimize & kAltMinos))
    optimize |= kMinosData;

  if (!doCond) optimizeTS= "U"+optimizeTS;

  if (fileName.empty()) {
    SetDefaults();
    Usage();
    return 1;
  }

  return 0;
}
