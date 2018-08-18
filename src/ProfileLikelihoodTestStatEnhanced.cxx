// @(#)root/roostats:$Id: ProfileLikelihoodTestStatEnhanced.cxx 754546 2016-06-13 18:07:15Z adye $
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additional Contributions: Giovanni Petrucciani
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#define USE_LASTHESSIAN

#include "Utilities/WorkspaceCalculatorConfig.h"
#include "ProfileLikelihoodTestStatEnhanced.h"
#include "RooFitResult.h"
#include "RooPullVar.h"
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,0)
#include "RooStats/DetailedOutputAggregator.h"
#define USE_DetailedOutputAggregator
#endif

#include "RooProfileLL.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
#include "TStopwatch.h"
#include "RooRandom.h"
#include "RooWorkspace.h"

#include "RooStats/HypoTestCalculatorGeneric.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ModelConfig.h"

#include <vector>
#include <limits>
#include "TH1.h"
#include "TTree.h"
#include "Math/DistFunc.h"
#ifdef USE_LASTHESSIAN
#include "TMethodCall.h"
#include "TMatrixDSym.h"
#endif

#include "WorkspaceCalculator.h"

#if (ROOT_VERSION_CODE==ROOT_VERSION(5,34,8) && ROOT_BUILD_SVN_REVISION >= 534080007) || \
     ROOT_VERSION_CODE>=ROOT_VERSION(5,34,9)
#define USE_EVALCOUNTER
#endif

using std::cout;
using std::cerr;
using std::endl;

#ifndef DEFINED_NAN
#define DEFINED_NAN
static Double_t MakeNaN() { return std::numeric_limits<Double_t>::quiet_NaN(); }
static const Double_t NaN = MakeNaN();
#endif

Bool_t RooStats::ProfileLikelihoodTestStatEnhanced::fgAlwaysReuseNll = kTRUE ;

void RooStats::ProfileLikelihoodTestStatEnhanced::SetAlwaysReuseNLL(Bool_t flag) { fgAlwaysReuseNll = flag ; }

const ULong64_t RooStats::ProfileLikelihoodTestStatEnhanced::default_optimize = (
   RooStats::ProfileLikelihoodTestStatEnhanced::kAdjustInitialMu    |
   RooStats::ProfileLikelihoodTestStatEnhanced::kResetNPOnNaN       |
   RooStats::ProfileLikelihoodTestStatEnhanced::kDoInitialScanOnNaN |
   RooStats::ProfileLikelihoodTestStatEnhanced::kSkipFitOnNaN       |
   RooStats::ProfileLikelihoodTestStatEnhanced::kRetryUncond        |
   RooStats::ProfileLikelihoodTestStatEnhanced::kZeroMu             |
   RooStats::ProfileLikelihoodTestStatEnhanced::kPrintResult        |
   RooStats::ProfileLikelihoodTestStatEnhanced::kConstOpt           |
   RooStats::ProfileLikelihoodTestStatEnhanced::kCheckMinimization  |
   RooStats::ProfileLikelihoodTestStatEnhanced::kReuseMinimizer     |
   RooStats::ProfileLikelihoodTestStatEnhanced::kCanRetryFit        |
   RooStats::ProfileLikelihoodTestStatEnhanced::kFallbackSimplex    |
   RooStats::ProfileLikelihoodTestStatEnhanced::kUseNllVal
);

struct RooStats::ProfileLikelihoodTestStatEnhanced::TreeVars {
  Double_t weight, ml, cml, pll, pllcmp, pll_init, mu_gen, mu_init, mu_init_err, fit_time, fit_time_cond, fit_time_retry;
  Double_t mu_fit[3], ml_fit[3], edm[3];
  Int_t is_data, is_alt, sample, have_uncond, have_cond, status, status_cond, status_init, tries, tries_cond, tries_retry;
  Int_t nfcn, nfcn_cond, nfcn_retry;
  Int_t covqual[3];
  std::vector<Double_t> uvars, uerrs, ivars, ierrs, cvars, cerrs, gvars;
  std::vector<std::string> names, gnames;
};

RooStats::ProfileLikelihoodTestStatEnhanced::ProfileLikelihoodTestStatEnhanced() {
       // Proof constructor. Do not use.
       fPdf = 0;
       fNll = fNll0 = 0;
       fCachedBestFitParams = 0;
       fLastData = 0;
       fLimitType = twoSided;
       fSigned = false;
       fDetailedOutputWithErrorsAndPulls = false;
       fDetailedOutputEnabled = 0;
       fDetailedOutput = NULL;

       fVarName = "Profile Likelihood Ratio";
       fReuseNll = false;
       fMinimizer=::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
       fStrategy=::ROOT::Math::MinimizerOptions::DefaultStrategy();
       fTolerance=TMath::Max(1.,::ROOT::Math::MinimizerOptions::DefaultTolerance());
       fTolerance0=-1.0;
       fPrintLevel=::ROOT::Math::MinimizerOptions::DefaultPrintLevel();

       fFillTree = false;
       fTree = 0;
       fVars = 0;
       fFlagData = 0;
       fSkipToys = 0;
       fComparisonTestStat = 0;
       fMinos = false;
       fMinosSet = 0;
       fInfo.SetOwner(kTRUE);
       fToyNum = fEvalNum = 0;
       fSeed = 0;
       fParmProf = -1;
       fSaveWS = 0;
       fSaveGlob = 0;
       fSavedInitialSnapshot = fUseSaved = false;
       fProblems = 0;
       fOptimize = default_optimize;
       fInitialValues = fInitialValuesCond = fFixInitially = 0;
       fCalc = 0;
       resultU = resultC = 0;
       fMinim = 0;
       fLastInitalSnap[0] = fLastInitalSnap[1] = fLastInitalSnap[2] = 0;
       fUseNll = false;
}

RooStats::ProfileLikelihoodTestStatEnhanced::ProfileLikelihoodTestStatEnhanced(RooAbsPdf& pdf) {
       fPdf = &pdf;
       fNll = fNll0 = 0;
       fCachedBestFitParams = 0;
       fLastData = 0;
       fLimitType = twoSided;
       fSigned = false;
       fDetailedOutputWithErrorsAndPulls = false;
       fDetailedOutputEnabled = 0;
       fDetailedOutput = NULL;

       fVarName = "Profile Likelihood Ratio";
       fReuseNll = false;
       fMinimizer=::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
       fStrategy=::ROOT::Math::MinimizerOptions::DefaultStrategy();
       // avoid default tolerance to be too small (1. is default in RooMinimizer)
       fTolerance=TMath::Max(1.,::ROOT::Math::MinimizerOptions::DefaultTolerance());
       fTolerance0=-1.0;
       fPrintLevel=::ROOT::Math::MinimizerOptions::DefaultPrintLevel();
       fFillTree = false;
       fTree = 0;
       fVars = 0;
       fFlagData = 0;
       fSkipToys = 0;
       fComparisonTestStat = 0;
       fMinos = false;
       fMinosSet = 0;
       fInfo.SetOwner(kTRUE);
       fToyNum = fEvalNum = 0;
       fSeed = 0;
       fParmProf = -1;
       fSaveWS = 0;
       fSaveGlob = 0;
       fSavedInitialSnapshot = fUseSaved = false;
       fProblems = 0;
       fOptimize = default_optimize;
       fInitialValues = fInitialValuesCond = fFixInitially = 0;
       fCalc = 0;
       resultU = resultC = 0;
       fMinim = 0;
       fLastInitalSnap[0] = fLastInitalSnap[1] = fLastInitalSnap[2] = 0;
       fUseNll = false;
}

RooStats::ProfileLikelihoodTestStatEnhanced::~ProfileLikelihoodTestStatEnhanced() {
       delete fMinim;
       if(fNll0!=fNll) delete fNll0;
       if(fNll) delete fNll;
       if(fCachedBestFitParams) delete fCachedBestFitParams;
       if(fDetailedOutput) delete fDetailedOutput;
       delete fTree;
       delete fVars;
       delete fComparisonTestStat;
       if (!(fOptimize & kSaveFitResult)) {
         delete resultU;
         delete resultC;
       }
       delete fLastInitalSnap[0];
       delete fLastInitalSnap[1];
       delete fLastInitalSnap[2];
}

TTree* RooStats::ProfileLikelihoodTestStatEnhanced::GetTree(Bool_t getOwnership)  {
  TTree* t=fTree;
  if (getOwnership && t) {
    t->ResetBranchAddresses();
    fTree=0;
  }
  return t;
}

ULong64_t RooStats::ProfileLikelihoodTestStatEnhanced::SetOptimize(const char* a, int verbose) {
  // optimize flag characters must correspond to optimize_names and OptimizeBits enums.
  static const char* optimize_chars = "DSrafNnIixuXzRpUCFBc10mOHPMvedsVw";
  static const char* optimize_names[]= {
    "SkipData",
    "SaveSeeds",
    "SaveFitResult",
    "AdjustInitialMu",
    "RetrySimplex",
    "ResetNP",
    "ResetNPOnNaN",
    "DoInitialScan",
    "DoInitialScanOnNaN",
    "SkipFitOnNaN",
    "RetryUncond",
    "NoFit",
    "ZeroMu",
    "PrintResult",
    "SaveInitialNP",
    "OnlyUnconditional",
    "OnlyConditional",
    "NoRestoreFit",
    "PruneBinning",
    "ConstOpt",
    "ConstOpt1",
    "DoOffset",
    "CheckMinimization",
    "SeparateOffsetNll",
    "Hesse",
    "PrintFittedParms",
    "ReuseMinimizer",
    "GetInitialNllVal",
    "CanRetryFit",
    "CloneData",
    "FallbackSimplex",
    "UseNllVal",
    "kConstConstraints",
  };

  if (!a) return fOptimize;
  bool neg= false;
  for (;*a; a++) {
    if      (*a=='-') neg= true;
    else if (*a=='+') neg= false;
    else if (const char* ac= strchr (optimize_chars, *a)) {
      if (neg) fOptimize &= ~(1LL<<(ac-optimize_chars));
      else     fOptimize |=  (1LL<<(ac-optimize_chars));
    } else {
      cerr << "ProfileLikelihoodTestStatEnhanced: invalid option: " << *a << endl; break;
    }
  }

  if (fOptimize & kConstOpt1) fOptimize &= ~kConstOpt;  // for clarity in debug message

  if (verbose<0 ? (fPrintLevel>0) : verbose) {
    cout << "ProfileLikelihoodTestStatEnhanced options: ";
    for (int i=0, j=0, l=strlen(optimize_chars); i<l; i++) {
      if (fOptimize & (1LL<<i)) cout << (j++ ? ", " : "") << optimize_names[i];
    }
    cout << endl;
  }
  return fOptimize;
}

TString RooStats::ProfileLikelihoodTestStatEnhanced::Join (const RooAbsCollection& c, const TString& sep, bool showname)
{
  // Borrow from WorkspaceCalculator
  return WorkspaceCalculator::Join(c,sep,showname);
}

RooArgSet* RooStats::ProfileLikelihoodTestStatEnhanced::GetCommonSnapshot (const RooAbsCollection& a, const RooAbsCollection& b)
{
  // Get a snapshot of elements of a that are in b, ordered as they are in b.
  // This is like RooAbsCollection::selectCommon + snapshot, but that orders according to a.
  RooArgSet s;
  for (RooFIter it= b.fwdIterator(); const RooAbsArg* vb= it.next();) {
    if (RooAbsArg* va= a.find(*vb)) s.add(*va);
  }
  return dynamic_cast<RooArgSet*>(s.snapshot());
}

TMatrixDSym RooStats::ProfileLikelihoodTestStatEnhanced::reducedCovarianceMatrix (const RooFitResult* res, const RooArgList& params)
{
  // Fixed version of RooFitResult::reducedCovarianceMatrix [see ROOT-8044] (without const params warning).
  // Return a reduced covariance matrix (Note that Vred _is_ a simple sub-matrix of V,
  // row/columns are ordered to matched the convention given in input argument 'params')

  const RooArgList& finalPars= res->floatParsFinal();

  std::vector<Int_t> ind;
  for (RooFIter it= params.fwdIterator(); const RooAbsArg* a= it.next();) {
    Int_t i= finalPars.index(a->GetName());
    if (i>=0) ind.push_back(i);
  }

  const TMatrixDSym& V= res->covarianceMatrix();
  Int_t n= ind.size();
  TMatrixDSym Vred(n);
  for (Int_t i=0; i<n; i++)
    for (Int_t j=0; j<n; j++)
      Vred(i,j)= V(ind[i], ind[j]);

  return Vred;
}

void RooStats::ProfileLikelihoodTestStatEnhanced::TCollectionRemoveAll (TCollection& c, const char* name)
{
  // Remove all objects with the given option or name
  bool cont;
  do {
    cont= false;
    for (TIter it= &c; TObject* o= it();) {
      const char* n= it.GetOption();
      if (!*n) n= o->GetName();
      if (strcmp (n, name) == 0) {
        c.Remove(o);
        cont= true;
        break;  // start again in case TIter was invalidated
      }
    }
  } while (cont);
}



Double_t RooStats::ProfileLikelihoodTestStatEnhanced::EvaluateProfileLikelihood(int type, RooAbsData& dataIn, RooArgSet& paramsOfInterest) {
        // interna function to evaluate test statistics
        // can do depending on type:
        // type  = 0 standard evaluation, type = 1 find only unconditional NLL minimum, type = 2 conditional MLL

  cout << "XXX start of EvaluateProfileLIkelihood " << endl;
       if( fDetailedOutputEnabled && fDetailedOutput ) {
	       delete fDetailedOutput;
	       fDetailedOutput = 0;
       }
       if( fDetailedOutputEnabled && !fDetailedOutput ) {
	       fDetailedOutput = new RooArgSet();
       }
       if (!(fOptimize & kSaveFitResult)) {
         delete resultU;
         delete resultC;
       }
       resultU= resultC= 0;
       if (!&dataIn) {
	 cout << "problem with data" << endl;
	 return 0 ;
       }

       fProblems = 0;

       RooFit::MsgLevel oldMsglevel = RooMsgService::instance().globalKillBelow();
       Int_t is_data=0, is_alt=-1;
       const RooArgSet* genPOIs= 0;
       Double_t mu_gen= NaN;
       if (fFlagData && &dataIn==fFlagData) is_data = 1;
       if (fCalc) {
         // work out what sort of fit we are doing from the HypoTestCalculatorGeneric
         // that we were told about.
         if (!fFlagData && &dataIn==fCalc->GetData()) is_data = 1;
         if (!is_data)
           if (RooStats::ToyMCSampler* toymcs= dynamic_cast<RooStats::ToyMCSampler*>(fCalc->GetTestStatSampler()))
             if (strcmp (fCalc->GetNullModel()->GetName(), fCalc->GetAlternateModel()->GetName()) != 0) {
               std::string mcname= toymcs->GetSamplingDistName();
               if      (mcname == fCalc->GetNullModel()     ->GetName()) is_alt=0;
               else if (mcname == fCalc->GetAlternateModel()->GetName()) is_alt=1;
               if (is_alt!=-1) {
                 const RooStats::ModelConfig* mc= is_alt ? fCalc->GetAlternateModel() : fCalc->GetNullModel();
                 genPOIs= mc->GetSnapshot();
                 if (const RooRealVar* mu_genVar= dynamic_cast<const RooRealVar*>(genPOIs->first()))
                   mu_gen = mu_genVar->getVal();
               }
             }
       }

       if ((fOptimize & kSaveSeeds) && (fToyNum==0 || !is_data)) {
         // save initial seed (which we see when we do the data) to toy_0, others to toy_1,...
         fInfo.Add (RooRandom::randomGenerator()->Clone(Form("rnd/toy_%d", is_data ? 0 : fToyNum+1)));
       }
       if (fToyNum==0 && is_data) fSeed= RooRandom::randomGenerator()->GetSeed();

       RooAbsData* readData= 0;
       if (fSaveWS && !is_data) {
         TString snap   = Form("toyGlob_%d",fToyNum);
         TString dsname = Form("toyData_%d",fToyNum);
         if (fUseSaved) {
           RooMsgService::instance().setGlobalKillBelow (RooFit::FATAL);
           Bool_t ok= fSaveWS->loadSnapshot(snap);
           RooMsgService::instance().setGlobalKillBelow(oldMsglevel);
           readData= fSaveWS->data(dsname);
           if (!readData && !ok) {
             cerr << "No toy";
             if (!readData) cerr << " dataset '"<<dsname<<"'";
             if (!readData && !ok) cerr << " or";
             if (!ok)       cerr << " global observable snapshot '"<<snap<<"'";
             cerr << " in workspace '"<<fSaveWS->GetName()<<"'"<<endl;
             return NaN;
           }
           if (fPrintLevel > 0) cout << "Read toy dataset '"<<dsname<<"' and global observable snapshot '"<<snap<<"' from workspace '"<<fSaveWS->GetName()<<"'"<<endl;
         } else {
           if (!fInfo.FindObject(fSaveWS)) fInfo.Add(fSaveWS);
           if (fSaveGlob) {
             fSaveWS->saveSnapshot(snap,*fSaveGlob);
           }
           RooAbsData* data2= dynamic_cast<RooAbsData*>(dataIn.Clone(dsname));
           Bool_t err= fSaveWS->import (*data2);
           delete data2;
           if (err)
             cerr << "Could not save toy dataset '"<<dsname<<"' to workspace '"<<fSaveWS->GetName()<<"'"<<endl;
           else if (fPrintLevel > 0) {
             cout << "Saved toy dataset '"<<dsname<<"'";
             if (fSaveGlob) cout << " and global observable snapshot '"<<snap<<"'";
             cout << " to workspace '"<<fSaveWS->GetName()<<"'"<<endl;
           }
         }
       }

       if (is_data) {
         if (fOptimize & kSkipData) return NaN;
       } else if (fSkipToys>0 && fToyNum<fSkipToys) {
         fToyNum++;
         return NaN;
       }

       int itype= 0;
       if (!((fOptimize & kOnlyUnconditional) && (fOptimize & kOnlyConditional))) {
         if (((fOptimize & kOnlyUnconditional) && type==2) ||
             ((fOptimize & kOnlyConditional)   && type==1)) {
           if (!is_data) fToyNum++;
           return NaN;
         }
         if (type==0) {
           if (fOptimize & kOnlyUnconditional) type= 1;
           if (fOptimize & kOnlyConditional)   type= 2;
           itype= type;
         }
       }

       int thisToy = fToyNum;
       if (is_data) fToyNum = -1;

       RooAbsData* prunedData= 0;
       if ((fOptimize & kPruneBinning) && !is_data) {
         prunedData= WorkspaceCalculator::SetDatasetBinning (fPdf, (readData?readData:&dataIn),
                                                             0, 0, 0, "ATLAS_H_*", fPrintLevel-2);
       }

       RooAbsData& data = prunedData ? *prunedData : readData ? *readData : dataIn;

       //data.Print("V");

       TStopwatch tsw;
       tsw.Start();

       RooArgSet* initialPOIs= dynamic_cast<RooArgSet*>(paramsOfInterest.snapshot());

       double initial_mu_value  = 0;
       RooRealVar* firstPOI = dynamic_cast<RooRealVar*>(paramsOfInterest.first());
       if (firstPOI) initial_mu_value = firstPOI->getVal();
       //paramsOfInterest.getRealValue(firstPOI->GetName());
       if (fPrintLevel > 1) {
            cout << "POIs: ";
            paramsOfInterest.Print("v");
       }

       RooFit::MsgLevel msglevel;
       if      (fPrintLevel < 2) msglevel= RooFit::FATAL;
       else if (fPrintLevel < 3) msglevel= RooFit::ERROR;
       else if (fPrintLevel < 4) msglevel= RooFit::INFO;
       else                      msglevel= RooFit::DEBUG;
       RooMsgService::instance().setGlobalKillBelow(msglevel);

       // make sure we set the variables that will be attached to this nll
       // This should be equivalent to fNll->getVariables(), but we can do it before we have the NLL
       if (fPrintLevel > 1) WorkspaceCalculator::PrintResourcesUsed("start getParameters");
       RooArgSet* attachedSet;
       if (fNll) attachedSet= fNll->getVariables();    // use old NLL for speed
       else      attachedSet= fPdf->getParameters(data);
       RooArgSet* allVars = dynamic_cast<RooArgSet*>(attachedSet->selectByAttrib("Constant",kFALSE));

       Double_t pllcmp=NaN;
       if (fComparisonTestStat) {
         const RooArgSet* saveAll = dynamic_cast<const RooArgSet*>(attachedSet->snapshot());
         RooMsgService::instance().setGlobalKillBelow(oldMsglevel);
         pllcmp = fComparisonTestStat->Evaluate(data,paramsOfInterest);
         RooMsgService::instance().setGlobalKillBelow(msglevel);
         *attachedSet = *saveAll;  // restore original parameters, which might have been changed by test stat Evaluate()
         delete saveAll;
       }

       // simple
       Bool_t reuse=(fReuseNll || fgAlwaysReuseNll) ;

       Bool_t created(kFALSE) ;
       if (!reuse || fNll==0) {
          delete fMinim; fMinim = 0;
          delete fNll;   fNll   = 0;
          RooArgSet* allParams = attachedSet;
          if (!(fOptimize & kConstConstraints)) {
            // Note that the NLL value will depend on which constraints are specified.
            // Since const parameters are excluded, the NLL value will jump if an NP is fixed or freed
            allParams = new RooArgSet(*allVars);
            allParams->add(paramsOfInterest,kTRUE);  // in case POIs are constraints and were set const
          }

          if (fPrintLevel > 1) {
            WorkspaceCalculator::PrintResourcesUsed();
            cout << "creating NLL from " << fPdf->ClassName() << "::" << fPdf->GetName() << "(" << fPdf
                 << ") and" << ((fOptimize & kCloneData) ? " cloned" : "") << " data "
                 << data.ClassName() << "::" << data.GetName() << "=" << data.numEntries() << "/" << data.sumEntries()
                 << " entries (" << &data << ") with " << allParams->getSize() << " constraints";
            if (fConditionalObs.getSize() > 0) cout << ", " << fConditionalObs.getSize() << " ConditionalObservables";
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,4)
            cout << ", with" << ((fOptimize & kDoOffset) ? "" : "out") << " offsetting";
#endif
            cout << endl;
          }

          // need to call constrain for RooSimultaneous until stripDisconnected problem fixed
          fNll = fPdf->createNLL(data, RooFit::CloneData(fOptimize & kCloneData),RooFit::Constrain(*allParams),RooFit::ConditionalObservables(fConditionalObs),
                                 RooFit::Verbose(fPrintLevel>2)
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,4)
                                 , RooFit::Offset(fOptimize & kDoOffset)
#endif
                                );
          if ((fOptimize & kSeparateOffsetNll))
            fNll0= fPdf->createNLL(data, RooFit::CloneData(fOptimize & kCloneData),RooFit::Constrain(*allParams),RooFit::ConditionalObservables(fConditionalObs),
                                   RooFit::Verbose(fPrintLevel>2));
          else
          fNll0= fNll;
          fUseNll= (fOptimize & kUseNllVal);

          created = kTRUE ;
          if (fPrintLevel > 1) {
            cout << "created  NLL " << fNll->ClassName() << "::" << fNll->GetName() << "(" << fNll;
            if (fNll0!=fNll) cout << " and " << fNll0;
            cout << ')' << endl;
            WorkspaceCalculator::PrintResourcesUsed();
          }
          if (!(fOptimize & kConstConstraints)) delete allParams;
       }
       if (reuse && !created) {
         if (&data != fLastData) {
           delete fMinim; fMinim = 0;
           if (fPrintLevel > 1) {
             cout << "reusing NLL " << fNll->ClassName() << "::" << fNll->GetName() << "(" << fNll
                  << ") with new" << ((fOptimize & kCloneData) ? " cloned" : "") << " data "
                  << data.ClassName() << "::" << data.GetName() << "=" << data.numEntries() << "/" << data.sumEntries()
                  << " (" << &data << ")" << endl;
           }
           fNll->setData (data, fOptimize & kCloneData) ;
           if (fNll0!=fNll) fNll0->setData (data, fOptimize & kCloneData) ;
           fUseNll= (fOptimize & kUseNllVal);
         } else if (fPrintLevel > 1) {
           cout << "reusing NLL " << fNll->ClassName() << "::" << fNll->GetName() << "(" << fNll << ")" << endl;
         }
       }


       if (fOptimize & kGetInitialNllVal) {
         Double_t nll0= fNll0->getVal();
         if (fPrintLevel > 0) cout << "EvaluateProfileLikelihood - Initial NLL = " << Form("%.4f",nll0) << endl;
         fUseNll= true;
       }

       int nVars = 0, ngVars = 0;
       if (fFillTree && !fTree) {
         fVars = new TreeVars;
         fTree = new TTree ("pll", "ProfileLikelihoodTestStatEnhanced statistics");
         fTree->Branch ("weight",      &fVars->weight);
         if (itype!=2) fTree->Branch ("ml",          &fVars->ml);
         if (itype!=1) fTree->Branch ("ml_cond",     &fVars->cml);
         fTree->Branch ("pll",         &fVars->pll);
         if (fComparisonTestStat) fTree->Branch ("pllcmp",  &fVars->pllcmp);
         if (fCalc||fFlagData)    fTree->Branch ("is_data", &fVars->is_data);
         if (itype==0) fTree->Branch ("pll_init",    &fVars->pll_init);
         if (fCalc)    fTree->Branch ("is_alt",      &fVars->is_alt);
         fTree->Branch ("sample",      &fVars->sample);
         if (fCalc)    fTree->Branch ("mu_gen",      &fVars->mu_gen);
         if (itype!=2) fTree->Branch ("have_uncond", &fVars->have_uncond);
         if (itype!=1) fTree->Branch ("have_cond",   &fVars->have_cond);
         if (itype!=2) fTree->Branch ("status",      &fVars->status);
         if (itype!=1) fTree->Branch ("status_cond", &fVars->status_cond);
         if (itype==0) fTree->Branch ("status_init", &fVars->status_init);
         if (itype!=2) fTree->Branch ("tries",       &fVars->tries);
         if (itype!=1) fTree->Branch ("tries_cond",  &fVars->tries_cond);
         if (itype==0) fTree->Branch ("tries_retry", &fVars->tries_retry);
         if (itype!=2) fTree->Branch ("fit_time",    &fVars->fit_time);
         if (itype!=1) fTree->Branch ("fit_time_cond",&fVars->fit_time_cond);
         if (itype==0) fTree->Branch ("fit_time_retry",&fVars->fit_time_retry);
#ifdef USE_EVALCOUNTER
         if (itype!=2) fTree->Branch ("nfcn",        &fVars->nfcn);
         if (itype!=1) fTree->Branch ("nfcn_cond",   &fVars->nfcn_cond);
         if (itype==0) fTree->Branch ("nfcn_retry",  &fVars->nfcn_retry);
#endif
         if (fOptimize & kGetFitResult) {
         if (itype!=2) fTree->Branch ("covqual",     &fVars->covqual[0]);
         if (itype!=1) fTree->Branch ("covqual_cond",&fVars->covqual[1]);
         if (itype==0) fTree->Branch ("covqual_init",&fVars->covqual[2]);
         if (itype!=2) fTree->Branch ("edm",         &fVars->edm[0]);
         if (itype!=1) fTree->Branch ("edm_cond",    &fVars->edm[1]);
         if (itype==0) fTree->Branch ("edm_init",    &fVars->edm[2]);
         }
         fTree->Branch ("toynum",      &fToyNum);
         fTree->Branch ("seed",        &fSeed);
         if (itype!=2) fTree->Branch ("mu_fit",      &fVars->mu_fit[0]);
         if (itype==0) fTree->Branch ("mu_fit_init", &fVars->mu_fit[2]);
         if (itype!=2) fTree->Branch ("ml_fit",      &fVars->ml_fit[0]);
         if (itype!=1) fTree->Branch ("ml_fit_cond", &fVars->ml_fit[1]);
         if (itype==0) fTree->Branch ("ml_fit_init", &fVars->ml_fit[2]);
         if (itype==0) fTree->Branch ("mu_init",     &fVars->mu_init);
         if (itype==0) fTree->Branch ("mu_init_err", &fVars->mu_init_err);
         fTree->SetDirectory(0);   // memory resident ntuple
         for (RooLinkedListIter it= paramsOfInterest.iterator(); TObject* vo= it.Next();) { // paramsOfInterest first
           RooRealVar* v= dynamic_cast<RooRealVar*>(vo);
           if (v) fVars->names.push_back(v->GetName());
         }
         for (RooLinkedListIter it= allVars->iterator(); TObject* vo= it.Next();) {
           RooRealVar* v= dynamic_cast<RooRealVar*>(vo);
           if (v && !paramsOfInterest.find(v->GetName())) fVars->names.push_back(v->GetName());
         }
         if (fSaveGlob) {
           for (RooLinkedListIter it= fSaveGlob->iterator(); TObject* vo= it.Next();) {
             if (RooRealVar* v= dynamic_cast<RooRealVar*>(vo)) fVars->gnames.push_back(v->GetName());
           }
         }
         nVars= fVars->names.size();
         ngVars= fVars->gnames.size();
         fVars->uvars.resize(nVars);
         fVars->uerrs.resize(nVars);
         fVars->cvars.resize(nVars);
         fVars->cerrs.resize(nVars);
         fVars->gvars.resize(ngVars);
         if (fOptimize & kSaveInitialNP) {
           fVars->ivars.resize(nVars);
           fVars->ierrs.resize(nVars);
         }
         if (itype!=2) {
           if (fPrintLevel > 1) cout << "EvaluateProfileLikelihood - NLL variables ("<<nVars<<"):";
           for (int i=0; i<nVars; i++) {
             if (fPrintLevel > 1) cout << ' ' << fVars->names[i];
             std::string ename= fVars->names[i]+"_err";
             fTree->Branch (fVars->names[i].c_str(), &fVars->uvars[i]);
             fTree->Branch (          ename.c_str(), &fVars->uerrs[i]);
           }
           if (fPrintLevel > 1) cout << endl;
         }
         if (fOptimize & kSaveInitialNP) {
           for (int i=0; i<nVars; i++) {
             std::string vname= fVars->names[i]+"_initial", ename= vname+"_err";
             fTree->Branch (vname.c_str(), &fVars->ivars[i]);
             fTree->Branch (ename.c_str(), &fVars->ierrs[i]);
           }
         }
         if (itype!=1) {
           for (int i=0; i<nVars; i++) {
             std::string vname= fVars->names[i]+"_cond", ename= vname+"_err";
             fTree->Branch (vname.c_str(), &fVars->cvars[i]);
             if (!paramsOfInterest.find(fVars->names[i].c_str())) fTree->Branch (ename.c_str(), &fVars->cerrs[i]);
           }
         }
         if (ngVars>0) {
           if (fPrintLevel > 1) cout << "EvaluateProfileLikelihood - Global Observables ("<<ngVars<<"):";
           for (int i=0; i<ngVars; i++) {
             if (fPrintLevel > 1) cout << ' ' << fVars->gnames[i];
             std::string gname= "globObs_"+fVars->gnames[i];
             fTree->Branch (gname.c_str(), &fVars->gvars[i]);
           }
           if (fPrintLevel > 1) cout << endl;
         }
       }
       if (fFillTree) {
         nVars= fVars->names.size();
         ngVars= fVars->gnames.size();
         fVars->have_uncond = fVars->have_cond = 0;
         for (int i=0; i<nVars; i++) {
           fVars->uvars[i] = fVars->cvars[i] = NaN;
           fVars->uerrs[i] = fVars->cerrs[i] = 0.0;
         }
         if (fOptimize & kSaveInitialNP) {
           for (int i=0; i<nVars; i++) {
             fVars->ivars[i] = NaN;
             fVars->ierrs[i] = 0.0;
           }
         }
         for (int i=0; i<3; i++) {
           fVars->mu_fit[i]= NaN;
           fVars->ml_fit[i]= NaN;
           fVars->covqual[i]= -3;
           fVars->edm[i]=    NaN;
         }
         fVars->weight= is_data ? 0.0 : 1.0;  // fill later from dataset
         fVars->is_data = is_data;
         fVars->is_alt  = is_alt;
         fVars->sample  = is_data ? -2 : 1;
         fVars->pllcmp  = pllcmp;
         fVars->mu_gen  = mu_gen;
         fVars->cvars[0]= initial_mu_value;   // in case we don't do the conditional fit
         for (int i=0; i<ngVars; i++) {
           if (RooRealVar* v= dynamic_cast<RooRealVar*>(fSaveGlob->find(fVars->gnames[i].c_str())))
             fVars->gvars[i] = v->getVal();
           else
             cerr << "EvaluateProfileLikelihood - Global Observable '"<<fVars->gnames[i]<<"' not found" << endl;
         }
       }

       if (!(fOptimize & kAdjustInitialMu))
         *attachedSet = paramsOfInterest;

       ///////////////////////////////////////////////////////////////////////
       // New profiling based on RooMinimizer (allows for Minuit2)
       // based on major speed increases seen by CMS for complex problems


       // other order
       // get the numerator
       RooArgSet* snap =  (RooArgSet*)paramsOfInterest.snapshot();

       if (fInitialValues && !is_data) {
         static int kilroy= 0;
         if (fPrintLevel>2 || (fPrintLevel>0 && !kilroy)) { cout << "Initialise fit with:" << endl; fInitialValues->Print("v"); }
         *allVars = *fInitialValues;
         allVars->setAttribAll("Constant",kFALSE);  // in case fInitialValues had some fixed vars
         kilroy++;
       }

       RooArgSet* origAttachedSet = (RooArgSet*) attachedSet->snapshot();

       if (fSaveWS && !is_data && !fUseSaved && !fSavedInitialSnapshot) {
         fSaveWS->saveSnapshot("initialParms",*attachedSet,true);
         fSavedInitialSnapshot = kTRUE;
       }

       tsw.Stop();
       double createTime = tsw.CpuTime();
       tsw.Start();

       // get the denominator
       double uncondML = 0, uncondML_init = 0;
       double fit_favored_mu = 0, muhat_init = 0, muhat_init_err = 0;
       int statusD = 0, triesD = 0, nfcnD = -1;
       int status_init = 0;
       RooArgSet* fittedPOIs = 0;
       if (type != 2) {
         RooRealVar* poi= dynamic_cast<RooRealVar*>(attachedSet->find(firstPOI->GetName()));
         if (poi) {
           if (!is_data && (fOptimize & kAdjustInitialMu) && !(fOptimize & kNoFit)) {
             // Choose a better starting position and error. Also MINUIT doesn't have to move mu back into range.
             Double_t mu_start = poi->getMin() + 1.05*poi->getError();
             if (poi->getVal() < mu_start && mu_start+1.05*poi->getError() < poi->getMax()) {
               if (fPrintLevel > 1)
                 cout << "Start unconditional fit at mu = " << mu_start << " +/- " << poi->getError()
                      << " (moved from " << poi->getVal() << ")" << endl;
               poi->setVal(mu_start);
             }
           }
         }

         if (fFillTree && (fOptimize & kSaveInitialNP)) {
           for (int i=0; i<nVars; i++) {
             if (RooRealVar* v= dynamic_cast<RooRealVar*>(attachedSet->find(fVars->names[i].c_str()))) {
               fVars->ivars[i] = v->getVal();
               fVars->ierrs[i] = v->hasError() && !v->isConstant() ? v->getError() : -1.0;
             } else
               cerr << "EvaluateProfileLikelihood - Initial variable '"<<fVars->names[i]<<"' not found" << endl;
           }
         }

         uncondML_init = uncondML = UnconditionalFit (1, statusD, triesD, nfcnD, attachedSet, &paramsOfInterest);
         status_init = statusD;
         // get best fit value for one-sided interval
         if (poi) {
           muhat_init = fit_favored_mu = poi->getVal();
           muhat_init_err = poi->hasError() ? poi->getError() : -1.0;
         }
         fittedPOIs= GetCommonSnapshot (*attachedSet, paramsOfInterest);
         if (is_data && !fInfo.FindObject("fitted_poi")) {
           RooAbsCollection* fitted_poi= fittedPOIs->snapshot();
           fitted_poi->setName("fitted_poi");
           fInfo.Add (fitted_poi);
         }

         if (fParmProf>=0 && fToyNum<fParmProf) {
           for (RooLinkedListIter it= allVars->iterator(); TObject* vo= it.Next();) {
             RooRealVar* v= dynamic_cast<RooRealVar*>(vo);
             if (!v) continue;
             TH1* h= fNll->createHistogram(v->GetName());
             h->SetDirectory(0);
             if (is_data) {
               h->SetName  (Form ("prof/%s_data",          v->GetName()));
               h->SetTitle (Form ("%s for data;%s;-ln(L)", v->GetName(), v->GetName()));
             } else {
               h->SetName  (Form ("prof/%s_toy%g_%d",               v->GetName(), mu_gen, fToyNum));
               h->SetTitle (Form ("%s for toy %d|#mu=%g;%s;-ln(L)", v->GetName(), fToyNum, mu_gen, v->GetName()));
             }
             fInfo.Add(h);
           }
         }
       }
       tsw.Stop();
       double fitTime1  = tsw.CpuTime();
       double fitTime3  = 0;

       //double ret = 0;
       int statusN = 0, triesN = 0, tries_retry = 0, nfcnN = -1, nfcn_retry = -1;
       tsw.Start();

       double condML = 0;

       bool doConditionalFit = (type != 1);

       // skip the conditional ML (the numerator) only when fit value is smaller than test value
       if (!fSigned && type==0 &&
           ((fLimitType==oneSided          && fit_favored_mu >= initial_mu_value) ||
            (fLimitType==oneSidedDiscovery && fit_favored_mu <= initial_mu_value))) {
          doConditionalFit = false;
          condML = uncondML;
          if (fFillTree) fVars->ml_fit[1]= fVars->ml_fit[0];
       }

       if ((fOptimize & kZeroMu) && type==0 && fit_favored_mu == initial_mu_value) {
          doConditionalFit = false;
          condML = uncondML;
          if (fFillTree) fVars->ml_fit[1]= fVars->ml_fit[0];
       }

       if (doConditionalFit) {


          RooArgSet* origNP= dynamic_cast<RooArgSet*>(origAttachedSet->selectByAttrib("Constant",kFALSE));
          origNP->remove(paramsOfInterest,kFALSE,kTRUE);

          //       cout <<" reestablish snapshot"<<endl;
          if (fOptimize & kResetNP) *attachedSet = *origNP;


          if (fInitialValuesCond) {
            *allVars = *fInitialValuesCond;
            allVars->setAttribAll("Constant",kFALSE);  // in case fInitialValuesCond had some fixed vars
          }

          *attachedSet = *snap;


          // set the POI to constant
          for (RooLinkedListIter it = paramsOfInterest.iterator(); TObject* tmpPar = it.Next();) {
            if (RooRealVar* tmpParA= dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName())))
              tmpParA->setConstant();
          }


          // in case no nuisance parameters are present
          // no need to minimize just evaluate the nll
          if (origNP->getSize() == 0) {
             fUseNll = true;
             condML = fNll0->getVal();
          }
          else {
             fNll->clearEvalErrorLog();
             RooFitResult* result = 0;
             condML = GetMinNLL(2, statusN, triesN, nfcnN, result, attachedSet, &paramsOfInterest,
                                (fOptimize & kResetNP) ? 0 : origNP);
#ifdef USE_DetailedOutputAggregator
             if( result && fDetailedOutputEnabled==1 ) {
               RooArgSet* detOutput = DetailedOutputAggregator::GetAsArgSet(result, "fitCond_", fDetailedOutputWithErrorsAndPulls);
               fDetailedOutput->addOwned(*detOutput);
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,3)
               delete detOutput;
#endif
             }
#endif
             if (!(fOptimize & kSaveFitResult)) delete resultC;
             resultC= result;

             if (fFillTree) {
               fVars->have_cond = 1;
               for (int i=0; i<nVars; i++) {
                 if (RooRealVar* v= dynamic_cast<RooRealVar*>(attachedSet->find(fVars->names[i].c_str()))) {
                   fVars->cvars[i] = v->getVal();
                   Double_t err = -1.0;
                   if (!v->isConstant()) {
                     if        (v->hasAsymError()) {
                       Double_t errlo= v->getAsymErrorLo(), errhi= v->getAsymErrorHi();
                       if      (errlo>=0.0) err =  errhi;
                       else if (errhi<=0.0) err = -errlo;
                       else                 err = 0.5 * (errhi - errlo);
                     } else if (v->hasError())
                       err = v->getError();
                   }
                   fVars->cerrs[i] = err;
                 } else
                   cerr << "EvaluateProfileLikelihood - Condiational variable '"<<fVars->names[i]<<"' not found" << endl;
               }
             }

             if ((fOptimize & kRetryUncond) && !(fOptimize & kNoFit) && type != 2 &&
                 ((statusD%1000!=0 && statusN%1000==0) || condML-uncondML < -0.025)) {
               fProblems++;
               tsw.Stop();
               TStopwatch tsw2;
               if (statusD%1000!=0 && statusN%1000==0)
                 cout << "    ----> conditional fit succeeded... will repeat failed unconditional fit" << endl;
               else
                 cout << "    ----> conditional fit is " << uncondML-condML << " better (ML "
                      << Form("%.4f",uncondML) << " at mu=" << fit_favored_mu
                      << " -> " << Form("%.4f",condML) << " at " << initial_mu_value
                      << ")... will repeat unconditional fit" << endl;
               for (RooLinkedListIter jt = paramsOfInterest.iterator(); TObject* tmpPar = jt.Next();) {
                 if (dynamic_cast<RooAbsArg*>(tmpPar)->isConstant()) continue;
                 if (RooRealVar* tmpParA= dynamic_cast<RooRealVar*>( attachedSet->find(tmpPar->GetName())))
                   tmpParA->setConstant(false);
               }
               int statusD2=0;
               double uncondML2 = UnconditionalFit (3, statusD2, tries_retry, nfcn_retry, attachedSet, &paramsOfInterest);
               double fit_favored_mu2 = 0.0;
               if (firstPOI) fit_favored_mu2 = attachedSet->getRealValue(firstPOI->GetName()) ;
               if (uncondML2>=uncondML) {
                 cout << "    ----> after re-fit, unconditional fit did not improve: ML "
                      << Form("%.4f",uncondML) << " at mu=" << fit_favored_mu <<  " -> "
                      << Form("%.4f",uncondML2) << " at " << fit_favored_mu2 << " - keep original" << endl;
               } else {
                 uncondML = uncondML2;
                 statusD = statusD2;
                 fit_favored_mu = fit_favored_mu2;
                 if (condML-uncondML < -0.025) {
                   cout << "    ----> after re-fit, unconditional fit is still worse than conditional: ML "
                        << Form("%.4f",uncondML) << " at mu=" << fit_favored_mu << " -> " << Form("%.4f",condML)
                        << " at " << initial_mu_value << endl;
                 }
                 if (!fSigned && type==0 &&
                     ((fLimitType==oneSided          && fit_favored_mu >= initial_mu_value) ||
                      (fLimitType==oneSidedDiscovery && fit_favored_mu <= initial_mu_value))) {
                   cout << "    ----> after re-fit, muhat " << fit_favored_mu
                        << (fLimitType==oneSided ? " > " : " < ") << initial_mu_value
                        << ", test statistic will be zero." << endl;
                   doConditionalFit = false;
                 }
               }
               tsw2.Stop();
               fitTime3 = tsw2.CpuTime();
               tsw.Start(kFALSE);
             }
          }
          delete origNP;

       }

       tsw.Stop();
       double fitTime2 = tsw.CpuTime();

       double pll, pll_init;
       if      (type == 1) pll_init = pll = uncondML;
       else if (type == 2) pll_init = pll = condML;
//     else if (statusD%1000 || statusN%1000)
//       pll_init = pll = -1.0;
       else if (!doConditionalFit)
         pll_init = pll = 0.0;
       else {
         pll      = condML-uncondML;
         pll_init = condML-uncondML_init;
         if (fSigned) {
           if (pll<0.0) pll = 0.0;   // bad fit
           if (fLimitType==oneSidedDiscovery ? (fit_favored_mu < initial_mu_value)
                                             : (fit_favored_mu > initial_mu_value))
             pll = -pll;
         }
       }

       if( fDetailedOutputEnabled ) {
         fDetailedOutput->addOwned(*new RooRealVar("ntentry","stats ntuple entry", fFillTree ? fTree->GetEntries() : -1));
         if (fPrintLevel > 2) {
           cout << "Detailed output:" << endl;
           fDetailedOutput->Print("v");
         }
       }

       if (fFillTree) {
         fVars->ml=  uncondML;
         fVars->cml= condML;
         fVars->pll= pll;
         fVars->pll_init= pll_init;
         fVars->mu_init= muhat_init;
         fVars->mu_init_err= muhat_init_err;
         fVars->status = statusD;
         fVars->status_cond = statusN;
         fVars->status_init = status_init;
         fVars->tries = triesD;
         fVars->tries_cond = triesN;
         fVars->tries_retry = tries_retry;
         fVars->fit_time = fitTime1;
         fVars->fit_time_cond = fitTime2;
         fVars->fit_time_retry = fitTime3;
         fVars->nfcn = nfcnD;
         fVars->nfcn_cond = nfcnN;
         fVars->nfcn_retry = nfcn_retry;
         fTree->Fill();
       }

       if (is_data && !fInfo.FindObject("fit_result") && !fInfo.FindObject("fit_result_cond")) {
         if (resultU) fInfo.Add(resultU->Clone("fit_result"));
         if (resultC) fInfo.Add(resultC->Clone("fit_result_cond"));
       }

       if (fPrintLevel > -1) {
          std::cout << "EvaluateProfileLikelihood - mu " << Join(*initialPOIs);
          if (is_data)
            std::cout << "|data";
          else if (genPOIs && genPOIs->getSize()>0)
            std::cout << "|" << Join(*genPOIs) << " toy " << fToyNum;
          if (type != 2)
            std::cout << " muhat " << Join(*fittedPOIs) <<  " ML " << Form("%.4f",uncondML);
          if (type != 1)
            std::cout << " cond-ML " << Form("%.4f",condML);
          if (type == 0)
            std::cout << " pll " << pll;
          if (fComparisonTestStat)
            std::cout << " cmp " << pllcmp;
          std::cout << " create,fit1,2-time " << createTime << "," << fitTime1 << "," << fitTime2;
          if (statusD%1000 || statusN%1000 ||
              (type==0 && doConditionalFit && !(fOptimize & kNoFit) && condML-uncondML < -0.025))
            std::cout << " - FAILED " << statusD << "," << statusN << "," << condML-uncondML;
          std::cout << std::endl;
       }

       // need to restore the values ?
       if (!(fOptimize & kNoRestoreFit)) *attachedSet = *origAttachedSet;

       delete fittedPOIs;
       delete initialPOIs;
       delete genPOIs;
       delete allVars;
       delete attachedSet;
       delete origAttachedSet;
       delete snap;
       delete prunedData;

       if (!reuse) {
         delete fMinim; fMinim= 0;
         if (fNll0!=fNll) delete fNll0;
	 delete fNll;
	 fNll = fNll0 = 0;
         fUseNll = false;
       }

       RooMsgService::instance().setGlobalKillBelow(oldMsglevel);

       // if(statusN!=0 || statusD!=0)
	 // ret= -1; // indicate failed fit [ WVE commented since not used ]

       fToyNum= thisToy;
       if (!is_data) fToyNum++;
       fEvalNum++;
       fLastData= &data;

  cout << "XXX end of EvaluateProfileLIkelihood " << endl;
       return pll;

     }

double RooStats::ProfileLikelihoodTestStatEnhanced::UnconditionalFit(int stage, int& status, int& tries, int& nfcn,
                                                                     RooArgSet* attachedSet, const RooArgSet* pois)
{
  status = 0;
  // minimize and count eval errors
  fNll->clearEvalErrorLog();
  RooFitResult* result = 0;
  double uncondML = GetMinNLL(stage, status, tries, nfcn, result, attachedSet, pois);
#ifdef USE_DetailedOutputAggregator
  // save this snapshot
  if( result && fDetailedOutputEnabled==1 ) {
    RooArgSet* detOutput = DetailedOutputAggregator::GetAsArgSet(result, "fitUncond_", fDetailedOutputWithErrorsAndPulls);
    fDetailedOutput->addOwned(*detOutput);
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,3)
    delete detOutput;
#endif
  }
#endif
  if (!(fOptimize & kSaveFitResult)) delete resultU;
  resultU= result;

  if (fFillTree) {
    fVars->have_uncond = 1;
    int nVars= fVars->names.size();
    for (int i=0; i<nVars; i++) {
      if (RooRealVar* v= dynamic_cast<RooRealVar*>(attachedSet->find(fVars->names[i].c_str()))) {
        fVars->uvars[i] = v->getVal();
        fVars->uerrs[i] = v->hasError() && !v->isConstant() ? v->getError() : -1.0;
      } else
        cerr << "EvaluateProfileLikelihood - Uncondiational variable '"<<fVars->names[i]<<"' not found" << endl;
    }
  }
  return uncondML;
}

void RooStats::ProfileLikelihoodTestStatEnhanced::PrintFitStatus (int stage, const char* msg, int status,
                                                                  const RooMinimizer&
#ifdef USE_EVALCOUNTER
                                                                  minim
#endif
                                                                  , const RooFitResult* res,
                                                                  Double_t initialNll, Double_t initialPOI, const RooRealVar* poi,
                                                                  bool verbose)
{
  if (fPrintLevel <= (verbose||fProblems ? 0 : 1)) return;
  Double_t Nll = (fUseNll ? fNll0->getVal() : NaN);
  cout << Form("    ----> %s: stage %d, status %d", msg, stage, status);
#ifdef USE_EVALCOUNTER
  cout << ", Nfcn " << minim.evalCounter();
#endif
  if (res) cout << Form(", cov. quality %d, edm %.7f", res->covQual(), res->edm());
  if      (fUseNll) cout << Form(", nll %.4f -> %.4f (%+g)", initialNll, Nll, Nll - initialNll);
  else if (res)     cout << Form(", nll %.4f", res->minNll()) << ((fOptimize & kDoOffset) ? "+offset" : "");
  if (poi) {
    if (poi->isConstant())
      cout << Form(", %s %g", poi->GetName(), initialPOI);
    else
      cout << Form(", %s %g -> %g (%+g)",
                   poi->GetName(), initialPOI, poi->getVal(), poi->getVal() - initialPOI);
  }
  cout << endl;
  if (fPrintLevel>1) WorkspaceCalculator::PrintResourcesUsed();
}


bool RooStats::ProfileLikelihoodTestStatEnhanced::saveHessian (TList& info, const TString& hesseName, int verbose)
{
  // Save Minuit2's last Hessian (if available) in info.
#ifdef USE_LASTHESSIAN
  // call ROOT::Minuit2::MnHesse::lastHessian() if it is available in this ROOT version
  static TMethodCall lastHessianCall("ROOT::Minuit2::MnHesse::lastHessian","");
  Long_t lastHessianRet=0;
  lastHessianCall.Execute(lastHessianRet);
  TMatrixDSym* lastHessian= (TMatrixDSym*)lastHessianRet;
  if (lastHessian) {
    if (verbose>=1) cout << "Save Minuit2 Hessian as " << hesseName << endl;
    info.Add (lastHessian->Clone(), hesseName);
    return true;
  }
#endif
  return false;
}


double RooStats::ProfileLikelihoodTestStatEnhanced::GetMinNLL (int stage, int& status, int& tries, int& nfcn,
                                                               RooFitResult*& result,
                                                               RooArgSet* attachedSet, const RooArgSet* pois,
                                                               const RooArgSet* reset)
{
  int tries0= 0, nfcn0= -1;
  if (fFixInitially) {
    RooArgSet* vars = dynamic_cast<RooArgSet*>(attachedSet->selectByAttrib("Constant",kFALSE));
    RooAbsCollection* fixVars= vars->selectCommon(*fFixInitially);
    Int_t nfix= fixVars->getSize();
    if (nfix>0 && nfix<vars->getSize()) {
      fixVars->setAttribAll ("Constant", 1);
      GetMinNLL1 (10+(stage%10), status, tries0, nfcn0, result, attachedSet, pois, reset);
      fixVars->setAttribAll ("Constant", 0);
      if (!(fOptimize & kSaveFitResult)) delete result;
    }
    delete fixVars;
    delete vars;
  }
  double nll= GetMinNLL1 (stage, status, tries, nfcn, result, attachedSet, pois, reset);
  if (tries0>1) tries += tries0-1;
  if (nfcn0 >0) nfcn  += nfcn0;
  return nll;
}

// Originally based on robustMinimize from
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HiggsAnalysis/CombinedLimit/src/ProfiledLikelihoodRatioTestStatExt.cc?revision=1.18&view=markup (6 Jan 2012)
double RooStats::ProfileLikelihoodTestStatEnhanced::GetMinNLL1 (int stage, int& status, int& tries, int& nfcn,
                                                                RooFitResult*& result,
                                                                RooArgSet* attachedSet, const RooArgSet* poisnap,
                                                                const RooArgSet* reset)
{
    status = -1;
    tries = 0;
    result = 0;
    nfcn = -1;

    if (fOptimize & kNoFit) {
      status = 0;
      fUseNll = true;
      return fNll0->getVal();
    }

    TString resname, restitle;
    if (fToyNum<0) {
      resname .Form ("fitres/data_stage%d",                      stage);
      restitle.Form ("fit result for data, stage %d",            stage);
    } else {
      resname .Form ("fitres/toy%d_stage%d",            fToyNum, stage);
      restitle.Form ("fit result for toy %d, stage %d", fToyNum, stage);
    }

    TString minimizer = fMinimizer;
    TString algorithm = ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo();
    if ((fOptimize & kFallbackSimplex) && algorithm == "Migrad") algorithm = "Minimize"; // prefer to use Minimize instead of Migrad
    Int_t strategy = fStrategy;
    Double_t tol = (fTolerance0>=0.0 ? fTolerance0 : fTolerance);

    if (!(fOptimize & kReuseMinimizer)) { delete fMinim; fMinim= 0; }
    if (fMinim) {
#ifdef USE_EVALCOUNTER
      fMinim->zeroEvalCount();
#endif
      if (fPrintLevel>1) cout << "reusing";
    } else {
      fMinim= new RooMinimizer(*fNll);
      if (fPrintLevel>1) cout << "created";
    }
    
    Int_t optimizeConst= (fOptimize & kConstOpt1) ? 1 : (fOptimize & kConstOpt) ? 2 : 0;
    if (fPrintLevel>1) {
      cout << " RooMinimizer(" << minimizer << "," << algorithm << ",stategy=" << strategy
           << ",eps=" << tol << ",optimizeConst=" << optimizeConst << ") "
           << fMinim << " for NLL " << fNll->ClassName() << "::" << fNll->GetName() << "(" << fNll << ")"
           << endl;
    }
    RooMinimizer& minim= *fMinim;

    //LM: RooMinimizer.setPrintLevel has +1 offset - so subtruct  here -1 + an extra -2
    //TJA: PrintLevel -1 and 1 most useful
    minim.setPrintLevel((fPrintLevel <= 1) ? -1 : fPrintLevel -1);
    if (fPrintLevel>1) minim.setProfile(kTRUE);
    if (fPrintLevel>3) minim.setVerbose(kTRUE);
    minim.setStrategy(strategy);
    minim.setEps(tol);
    // this causes a memory leak
    minim.optimizeConst(optimizeConst);

#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,10)
    // if we set different default call limits in the MinimizerOptions, use those instead of RooMinimizer's default of 500*N.
    if (ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls() != 0)
      minim.setMaxFunctionCalls (ROOT::Math::MinimizerOptions::DefaultMaxFunctionCalls());
    if (ROOT::Math::MinimizerOptions::DefaultMaxIterations()    != 0)
      minim.setMaxIterations    (ROOT::Math::MinimizerOptions::DefaultMaxIterations());  // actually, MaxIterations isn't used by Minuit(2)
#endif

    double initialNll = (fUseNll ? fNll0->getVal() : NaN);
    double initialPOI = 0.0, poimin = 0.0;
    RooArgSet pois;
    for (RooFIter it= poisnap->fwdIterator(); const RooAbsArg* v= it.next();) {
      if (RooAbsArg* va= attachedSet->find(*v)) pois.add(*va);
    }
    RooRealVar* poivar= dynamic_cast<RooRealVar*>(pois.first());
    if (poivar) {
      initialPOI = poivar->getVal();
      poimin = poivar->getMin();
    }
    double lastNll= initialNll, lastPOI= initialPOI;

    RooArgSet* printVars = dynamic_cast<RooArgSet*>(attachedSet->selectByAttrib("Constant",kFALSE));
    RooArgSet cmpVars = *printVars;
    printVars->add (pois, kTRUE);  // add POIs back in in case they were constant (ie. conditional fit)
    printVars->sort();

    const char* stageName= (stage==1) ? "unconditional" :
                           (stage==2) ? "conditional"   :
                           (stage==3) ? "retried unconditional" : "";
    static int lastEvalNum= -1;
    if (stage==1 || fEvalNum!=lastEvalNum) {
      TString nllmsg= (fUseNll ? Form(" (NLL=%g)", initialNll) : "");
      if (fPrintLevel > 2) {
        cout << "EvaluateProfileLikelihood - Initial parameters before " << stageName << " fit" << nllmsg << ":" << endl;
        attachedSet->Print("v");
        cout << "EvaluateProfileLikelihood - end of initial parameters before " << stageName << " fit" << endl;
      } else if (fPrintLevel > 1 || (!fEvalNum && fPrintLevel > 0)) {
        int show=1;
        if (fPrintLevel > 1 && stage>=1 && stage<=3) {
          RooArgSet* snap= dynamic_cast<RooArgSet*>(cmpVars.snapshot());
          if (snap && fLastInitalSnap[stage-1]) {
            if (!fLastInitalSnap[stage-1]->equals(*snap))
              show= 2;
            else {
              show= 0;
              for (RooFIter it= snap->fwdIterator(); const RooAbsArg* a= it.next();) {
                const RooRealVar* v= dynamic_cast<const RooRealVar*>(a);
                const RooRealVar* l= dynamic_cast<const RooRealVar*>(fLastInitalSnap[stage-1]->find(*a));
                if (!v && !l) continue;
                if (!v || !l ||
                    v->getVal()   != l->getVal() ||
                    v->getError() != l->getError() ||
                    v->getMin()   != l->getMin() ||
                    v->getMax()   != l->getMax()) {
                  show= 2;
                  break;
                }
              }
            }
          }
          delete fLastInitalSnap[stage-1];
          fLastInitalSnap[stage-1]= snap;
        }
        if (show==0) {
          cout << "EvaluateProfileLikelihood - Same initial fit parameters before " << stageName << " fit" << nllmsg << endl;
        } else if (show==1) {
          cout << "EvaluateProfileLikelihood - Initial parameters before " << stageName << " fit" << nllmsg << ":" << endl;
          printVars->Print("v");
        } else if (show==2) {
          cout << "EvaluateProfileLikelihood - New initial parameters before " << stageName << " fit" << nllmsg << ":" << endl;
          printVars->Print("v");
        }
      }
    }
    lastEvalNum= fEvalNum;

    bool doInitialScan = (fOptimize & kDoInitialScan);
    if (fUseNll && (std::isnan(lastNll) || std::isinf(lastNll))) {
      fProblems++;
      bool stillBad = true;
      cout << "    ----> Stage " << stage << " initial NLL is " << lastNll << " at mu " << lastPOI;
      if ((fOptimize & kResetNPOnNaN)       && reset) {
        *attachedSet = *reset;
        initialNll = lastNll = fNll0->getVal();
        if (poivar) initialPOI = lastPOI = poivar->getVal();
        cout << " - resetting NPs gives NLL " << lastNll;
        stillBad = (std::isnan(lastNll) || std::isinf(lastNll));
      }
      if ((fOptimize & kDoInitialScanOnNaN) && stillBad) {
        cout << " - do initial Scan";
        doInitialScan = true;
        stillBad = false;
      }
      if ((fOptimize & kSkipFitOnNaN)       && stillBad) {
        cout << " - bail now!" << endl;
        status = -2;
        delete printVars;
        return lastNll;
      }
      cout << endl;
    }

    if (doInitialScan) {
      if (fPrintLevel>(fProblems?0:1)) cout << "    ----> initial Scan" << endl;
      status = minim.minimize(fMinimizer,"Scan");
      RooFitResult* res = 0;
      if (fOptimize & kGetFitResult)
        res = minim.save (Form("%s_try%d_scan",resname.Data(),tries),
                          Form("%s %s %s, try %d",fMinimizer.Data(),"Scan",restitle.Data(),tries));
      PrintFitStatus (stage, "scan result", status, minim, res, lastNll, lastPOI, poivar, false);
      if (fOptimize & kSaveFitResult) fInfo.Add(res);
      else                            delete res;
      if (fUseNll) initialNll = lastNll = fNll0->getVal();
      if (poivar) initialPOI = lastPOI = poivar->getVal();
      if (fUseNll && (fOptimize & kSkipFitOnNaN) && (std::isnan(lastNll) || std::isinf(lastNll))) {
        cout << "    ----> Stage " << stage << " after initial Scan, NLL is still " << lastNll << " at mu " << lastPOI
             << " - bail now!" << endl;
        status = -2;
#ifdef USE_EVALCOUNTER
        nfcn= minim.evalCounter();
#endif
        delete printVars;
        return lastNll;
      }
    }

    const int try1     = (fTolerance0>=0.0 ? 0 : 1);
    const int maxtries = 5-strategy;
    bool fillNllLater= false;
    if (fPrintLevel>1) WorkspaceCalculator::PrintResourcesUsed();
    if (fPrintLevel>(fProblems?0:1))
      cout << "    ----> Stage " << stage << ": " << algorithm << " with strategy " << strategy << " and tolerance " << tol << endl;
    for (tries = try1;; ++tries) {
        status = minim.minimize(minimizer, algorithm);
        if (status%1000 != 0) fProblems++;
        if (fDetailedOutputEnabled==1 || (fOptimize & kGetFitResult)) {
          if (!(fOptimize & kSaveFitResult)) delete result;
          result = minim.save (Form("%s_try%d",resname.Data(),tries),
                               Form("%s %s %s, try %d",minimizer.Data(),algorithm.Data(),restitle.Data(),tries));
          if (fOptimize & kSaveFitResult) {
            fInfo.Add(result);
            saveHessian (fInfo, Form("%s_try%d_hesse",resname.Data(),tries), fPrintLevel-1);
          }
          if (fFillTree && stage>=1 && stage<=3) {
            fVars->covqual[stage-1]= result->covQual();
            fVars->edm    [stage-1]= result->edm();
            if (stage==1) {
              fVars->covqual[2]= fVars->covqual[0];
              fVars->edm    [2]= fVars->edm    [0];
            }
          }
        }

        Double_t thisPOI = (poivar ? poivar->getVal() : NaN);
        Double_t thisNll = (fUseNll ? fNll0->getVal() : NaN);
        if (fFillTree && stage>=1 && stage<=3) {
          if (!fUseNll) fillNllLater= true;
          fVars->mu_fit[stage-1]= thisPOI;
          fVars->ml_fit[stage-1]= thisNll;
          if (stage==1) {
            fVars->mu_fit[2]= thisPOI;
            fVars->ml_fit[2]= thisNll;
          }
        }
        if ((fOptimize & kZeroMu) && (status%1000)==0 && poivar && !poivar->isConstant()) {
          Double_t poi0 = thisPOI-poimin;
          if (poi0>0.0 && poi0<1e-3 &&
              poivar->getErrorLo() < 0.0 && poi0 < -0.1*poivar->getErrorLo()) {
            poivar->setVal(poimin);
            Double_t newNll= (fUseNll ? fNll0->getVal() : NaN);
            if (fPrintLevel>(fProblems?0:1)) {
              cout << "    ----> Stage " << stage << ": moving mu " << thisPOI << " -> " << poimin;
              if (fUseNll) cout << " changed NLL " << Form("%.4f -> %.4f (%+g)",thisNll,newNll,newNll-thisNll);
            }
            if (!fUseNll || newNll<=thisNll) {
              if (fPrintLevel>(fProblems?0:1)) cout << " - use new setting" << endl;
            } else {
              if (fPrintLevel>(fProblems?0:1)) cout << " - keep old setting" << endl;
              poivar->setVal(thisPOI);
              if (fNll0->getVal() != thisNll) {
                cerr << "WARNING: restoring mu " << poimin << " -> " << thisPOI << " did not restore NLL: "
                     << Form("%.4f -> %.4f -> %.4f",thisNll,newNll,fNll0->getVal()) << endl;
              }
            }
          }
        }

        if        ((fOptimize & kCheckMinimization) && fUseNll &&
                   (status%1000) == 0 && fNll0->getVal() > initialNll + 0.02) {
            fProblems++;
            PrintFitStatus (stage, "false minimum", status, minim, result, lastNll, lastPOI, poivar);
            if (result) *attachedSet = result->floatParsInit();
            status = 1;
        } else if ((fOptimize & kCheckMinimization) &&
                   (status%1000) == 0 && result && result->covQual() == -1) {
            fProblems++;
            PrintFitStatus (stage, "cov fail", status, minim, result, lastNll, lastPOI, poivar);
            *attachedSet = result->floatParsInit();
            status = 1;
        } else if ((status%1000) == 0) {
            PrintFitStatus (stage, "success",       status, minim, result, lastNll, lastPOI, poivar, false);
        } else if (tries >= 2 && result && result->edm() < 0.05*tol &&
                   (!(fOptimize & kCheckMinimization) || !(fUseNll && fNll0->getVal() > initialNll + 0.02))) {
            PrintFitStatus (stage, "acceptable",    status, minim, result, lastNll, lastPOI, poivar);
            status = 0;
        } else if (tries < maxtries) {
            PrintFitStatus (stage, "partial fail",  status, minim, result, lastNll, lastPOI, poivar);
        } else {
            PrintFitStatus (stage, "final fail",    status, minim, result, lastNll, lastPOI, poivar);
        }
        if (tries>0 && ((status%1000) == 0 || !(fOptimize & kCanRetryFit) || tries >= maxtries)) break;

        if (fUseNll) lastNll = fNll0->getVal();
        if (poivar) lastPOI= poivar->getVal();

        if (tries == 0) {
          if (fOptimize & kDoOffset) {
            minim.setOffsetting(kFALSE);  // reset offset
            minim.setOffsetting(kTRUE);
          }
          tol= fTolerance;
          minim.setEps(tol);
        } else {
          // cout << "    ----> Doing a Scan before re-trying fit" << endl;
          status = minim.minimize(fMinimizer,"Scan");
          RooFitResult* res = 0;
          if (fOptimize & kGetFitResult)
            res = minim.save (Form("%s_try%d_scan",resname.Data(),tries),
                              Form("%s %s %s, try %d",fMinimizer.Data(),"Scan",restitle.Data(),tries));
          PrintFitStatus (stage, "scan result", status, minim, res, lastNll, lastPOI, poivar);
          if (fOptimize & kSaveFitResult) fInfo.Add(res);
          else                            delete res;
          if (fUseNll) lastNll = fNll0->getVal();
          if (poivar) lastPOI= poivar->getVal();
          if (tries+1 >= maxtries) {
            minimizer = "Minuit";
            algorithm = "migradimproved";
            strategy= 1;
            minim.setStrategy(strategy);
          } else if (tries > 1) {
            minim.setStrategy(++strategy);
          }
        }
        cout << "    ----> Re-trying " << algorithm << " with strategy " << strategy << " and tolerance " << tol << endl;
    }

    // Minimize already tries Simplex if Migrad fails, so don't do it again.
    if ((fOptimize & kRetrySimplex) && (status%1000) != 0) {
        if (fUseNll) lastNll = fNll0->getVal();
        if (poivar) lastPOI= poivar->getVal();
        cout << "    ----> Last attempt: simplex method" << endl;
        status = minim.minimize(fMinimizer,"Simplex");
        if (fDetailedOutputEnabled==1 || (fOptimize & kGetFitResult)) {
          if (!(fOptimize & kSaveFitResult)) delete result;
          result = minim.save (Form("%s_try%d_simplex",resname.Data(),tries),
                               Form("%s %s %s, try %d",fMinimizer.Data(),"Simplex",restitle.Data(),tries));
          if (fOptimize & kSaveFitResult) fInfo.Add(result);
        }
        if (!fUseNll || !(fNll0->getVal() > initialNll + 0.02)) {
            PrintFitStatus (stage, "success", status, minim, result, lastNll, lastPOI, poivar);
            status = 0;
        } else {
            PrintFitStatus (stage, "final fail", status, minim, result, lastNll, lastPOI, poivar);
            if ((status%1000) == 0) status = -1;
        }
    }

    if (fOptimize & kHesse) {
      if (fUseNll) lastNll = fNll0->getVal();
      if (poivar) lastPOI= poivar->getVal();
      if (fPrintLevel>1) cout << "    ----> Run HESSE" << endl;
      Int_t stat= minim.hesse();
      RooFitResult* res= 0;
      if (fDetailedOutputEnabled==1 || (fOptimize & kGetFitResult))
        res = minim.save (Form("%s_hesse",resname.Data()),
                          Form("%s %s %s",fMinimizer.Data(),"HESSE",restitle.Data()));
      PrintFitStatus (stage, "HESSE", stat, minim, result, lastNll, lastPOI, poivar, false);
      if (res) {
        if (fToyNum<0 && stage!=2 && pois.getSize()>=2) {
          TMatrixDSym* poicov= new TMatrixDSym (reducedCovarianceMatrix (res, pois));
          if (poicov->GetNrows()>=2) {  // could be less if some pois are constant
            TCollectionRemoveAll(fInfo,"fitted_poicov");
            fInfo.Add(poicov,"fitted_poicov");
          } else
            delete poicov;
        }
        if (fOptimize & kSaveFitResult) {
          fInfo.Add(res);
          if (stat==0) result= res;
          saveHessian (fInfo, Form("%s_hesse_hesse",resname.Data()), fPrintLevel-1);
        } else if (stat==0) {
          delete result;
          result= res;
        } else
          delete res;
      }
    }

    if (fMinos && fToyNum<0 &&
        stage == ((fOptimize & kOnlyConditional) ? 2 : 1)) {
      RooArgSet* vars = dynamic_cast<RooArgSet*>(attachedSet->selectByAttrib("Constant",kFALSE));
      Int_t stat= -999;
      RooArgSet mvars;
      RooArgSet& minos_results= fMinosSet ? mvars : *vars;
      if (fMinosSet) {
        RooArgSet xvars;
        for (RooFIter it= fMinosSet->fwdIterator(); const RooAbsArg* vm= it.next();) {
          if (RooAbsArg* v= vars->find(*vm)) mvars.add(*v);
          else                               xvars.add(*vm);
        }
        if (xvars.getSize() > 0) {
          cout << "MINOS parameters not defined: ";
          xvars.Print();
        }
        if (mvars.getSize() > 0) {
          if (fPrintLevel>0) {
            cout << "    ----> Run MINOS for ";
            mvars.Print();
          }
          stat= minim.minos(mvars);
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,34,10) && ROOT_VERSION_CODE<=ROOT_VERSION(5,34,14)
          minim.fitter()->Config().SetMinosErrors(false);  // don't rerun MINOS if RooMinimizer is reused
#endif
        } else
          cout << "    ----> MINOS: nothing to do";
      } else {
        if (fPrintLevel>0) cout << "    ----> Run MINOS" << endl;
        vars->Print();
        stat= minim.minos();
      }
      if (stat != -999) {
        if (fPrintLevel > 0) {
          cout << "MINOS status " << stat << ", results:-" << endl;
          minos_results.Print("v");
        } else if (stat)
          cout << "MINOS status " << stat << endl;
        if (!fInfo.FindObject("minos_results")) {
          RooAbsCollection* snap= minos_results.snapshot();
          snap->setName("minos_results");
          fInfo.Add(snap);
        }
        if (fOptimize & kSaveFitResult) {
          RooFitResult* res = minim.save (Form("%s_minos",resname.Data()),
                                          Form("%s %s %s",fMinimizer.Data(),"MINOS",restitle.Data()));
          fInfo.Add(res);
        }
        if (fPrintLevel>1) WorkspaceCalculator::PrintResourcesUsed();
      }
      delete vars;
    }

    lastNll = fNll0->getVal();
    fUseNll = true;

    if (fillNllLater) {
      fVars->ml_fit[stage-1]= lastNll;
      if (stage==1) fVars->ml_fit[2]= lastNll;
    }

    if ((status%1000) == 0 && (std::isnan(lastNll) || std::isinf(lastNll))) {
      if (poivar) lastPOI= poivar->getVal();
      cout << "    ----> Stage " << stage << " fitted NLL is " << lastNll << " at mu " << lastPOI << endl;
      status= -2;
    }
#ifdef USE_EVALCOUNTER
    nfcn= minim.evalCounter();
#endif

    if ((fOptimize & kPrintFittedParms) || fPrintLevel > 2 || (!fEvalNum && fPrintLevel > 1)) {
      cout << "EvaluateProfileLikelihood - Parameters after " << stageName << " fit (NLL=" << lastNll << "):" << endl;
      printVars->Print("v");
    } else if (stage==1 && fToyNum<0 && fPrintLevel > 0) {
      cout << "EvaluateProfileLikelihood - POI after " << stageName << " fit to data (NLL=" << lastNll << "):" << endl;
      pois.Print("v");
    }
    delete printVars;

    return lastNll;
}
