// @(#)root/roostats:$Id: ProfileLikelihoodTestStatEnhanced.h 754546 2016-06-13 18:07:15Z adye $
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Additional Contributions: Giovanni Petrucciani 
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ProfileLikelihoodTestStatEnhanced
#define ROOSTATS_ProfileLikelihoodTestStatEnhanced

//_________________________________________________
/*
BEGIN_HTML
<p>
ProfileLikelihoodTestStatEnhanced is an implementation of the TestStatistic interface that calculates the profile
likelihood ratio at a particular parameter point given a dataset.  It does not constitute a statistical test, for that one may either use:
<ul>
 <li> the ProfileLikelihoodCalculator that relies on asymptotic properties of the Profile Likelihood Ratio</li>
 <li> the Neyman Construction classes with this class as a test statistic</li>
 <li> the Hybrid Calculator class with this class as a test statistic</li>
</ul>

</p>
END_HTML
*/
//

#include "TList.h"
#include "TMatrixDSym.h"
#include "RooArgSet.h"

#ifndef ROOSTATS_TestStatistic
#include "TestStatistic.h"
#endif

class TTree;
class RooMinimizer;
class RooAbsData;
class RooAbsPdf ;;
class RooAbsReal;
class RooRealVar;
class RooWorkspace;
class RooFitResult;

namespace RooStats {
  class HypoTestCalculatorGeneric;

  class ProfileLikelihoodTestStatEnhanced : public TestStatistic {

   private:
     struct TreeVars;
     enum LimitType {twoSided, oneSided, oneSidedDiscovery};

   public:

     // Fit optimization options. These must match the strings defined in SetOptimize().
     enum OptimizeBits {
       kSkipData           = 1LL <<  0, // D: skip calculation for data
       kSaveSeeds          = 1LL <<  1, // S: save all random number seeds for retrieval by GetInfo()
       kSaveFitResult      = 1LL <<  2, // r: save all fit results for retrieval by GetInfo()
       kAdjustInitialMu    = 1LL <<  3, // a: adjust mu if it is out of range before the unconditional fit
       kRetrySimplex       = 1LL <<  4, // f: Minimize already runs Simplex if Migrad fails, so dont do it again - unless enabled here
       kResetNP            = 1LL <<  5, // N: Reset NPs before conditional fit
       kResetNPOnNaN       = 1LL <<  6, // n: Reset NPs before conditional fit if the initial NLL is NaN or inf
       kDoInitialScan      = 1LL <<  7, // I: Do a Scan before trying Migrad
       kDoInitialScanOnNaN = 1LL <<  8, // i: Do a Scan before trying Migrad if the initial NLL is NaN or inf
       kSkipFitOnNaN       = 1LL <<  9, // x: Don't try Migrad if the initial NLL is NaN or inf
       kRetryUncond        = 1LL << 10, // u: Retry unconditional fit if TS<-0.025
       kNoFit              = 1LL << 11, // X: evaluate NLLs like SimpleLikelihoodRatioTestStat with alt=1 (may want to load snapshot)
       kZeroMu             = 1LL << 12, // z: allow mu=0 (or min val) as fit result
       kPrintResult        = 1LL << 13, // R: print fit results
       kSaveInitialNP      = 1LL << 14, // p: save nuisance parameters to ntuple before unconditional fit
       kOnlyUnconditional  = 1LL << 15, // U: skip   conditional fit
       kOnlyConditional    = 1LL << 16, // C: skip unconditional fit
       kNoRestoreFit       = 1LL << 17, // F: don't restore NP parameters after fit
       kPruneBinning       = 1LL << 18, // B: prune binning on llll dataset
       kConstOpt           = 1LL << 19, // c: enable constant term and tracking optimization (flag=2)
       kConstOpt1          = 1LL << 20, // 1: enable constant term and tracking optimization (flag=1)
       kDoOffset           = 1LL << 21, // 0: enable NLL offsetting
       kCheckMinimization  = 1LL << 22, // m: check that minimization didn't increase NLL
       kSeparateOffsetNll  = 1LL << 23, // O: use separate RooNLLVar without offsetting (with kDoOffset) as workaround for RooFit bug, fixed in 5.34 r49074. Safer to disable ReuseNLL (-O-R).
       kHesse              = 1LL << 24, // H: run Hesse after fit
       kPrintFittedParms   = 1LL << 25, // P: print fitted parameters after fit (same as with PrintLevel=2)
       kReuseMinimizer     = 1LL << 26, // M: reuse RooMinimizer object
       kGetInitialNllVal   = 1LL << 27, // v: get initial NLL value before creating RooMinimizer
       kCanRetryFit        = 1LL << 28, // e: if fit fails, scan and retry fit with different strategies (default)
       kCloneData          = 1LL << 29, // d: clone dataset when creating/updating RooNLLVar
       kFallbackSimplex    = 1LL << 30, // s: Fall back to Simplex if Migrad fails (ie. use "Minimize" instead of "Migrad")
       kUseNllVal          = 1LL << 31, // V: use NLL value before end
       kConstConstraints   = 1LL << 32, // w: apply constraints also to constant parameters
       kGetFitResult       = (kSaveFitResult | kPrintResult),
     };

     ProfileLikelihoodTestStatEnhanced();
     ProfileLikelihoodTestStatEnhanced(RooAbsPdf& pdf);
     virtual ~ProfileLikelihoodTestStatEnhanced();

     //LM use default copy constructor and assignment copying the pointers. Is this what we want ?

     void SetOneSided(Bool_t flag=true) {fLimitType = (flag ? oneSided : twoSided);}
     void SetOneSidedDiscovery(Bool_t flag=true) {fLimitType = (flag ? oneSidedDiscovery : twoSided);}
     void SetSigned(Bool_t flag=true) {fSigned = flag;}  // +/- t_mu instead of t_mu>0 with one-sided settings
     //void SetOneSidedDiscovery(Bool_t flag=true) {fOneSidedDiscovery = flag;}

     bool IsTwoSided() const { return fLimitType == twoSided; }
     bool IsOneSidedDiscovery() const { return fLimitType == oneSidedDiscovery; }

     static void SetAlwaysReuseNLL(Bool_t flag);

     void SetReuseNLL(Bool_t flag) { fReuseNll = flag ; }
     void SetLOffset(Bool_t flag=kTRUE) { if (flag) fOptimize |= kDoOffset; else fOptimize &= ~kDoOffset; }

     void SetMinimizer(const char* minimizer){ fMinimizer=minimizer;}
     void SetStrategy(Int_t strategy){fStrategy=strategy;}
     void SetTolerance(double tol){fTolerance=tol;}
     void SetTolerance0(double tol){fTolerance0=tol;}
     void SetPrintLevel(Int_t printlevel){fPrintLevel=printlevel;}
    
     // Main interface to evaluate the test statistic on a dataset
     virtual Double_t Evaluate(RooAbsData& data, RooArgSet& paramsOfInterest) {
        return EvaluateProfileLikelihood(0, data, paramsOfInterest);
     }

     // evaluate  the profile likelihood ratio (type = 0) or the minimum of likelihood (type=1) or the conditional LL (type = 2) 
     virtual Double_t EvaluateProfileLikelihood(int type, RooAbsData &data, RooArgSet & paramsOfInterest);
     
     virtual void EnableDetailedOutput( int e=1, bool withErrorsAndPulls=false ) {
        fDetailedOutputEnabled = e;
        fDetailedOutputWithErrorsAndPulls = withErrorsAndPulls;
        delete fDetailedOutput;
        fDetailedOutput = NULL;
     }
     virtual const RooArgSet* GetDetailedOutput(void) const {
	     // Returns detailed output. The value returned by this function is updated after each call to Evaluate().
	     // The returned RooArgSet contains the following:
	     // <ul>
	     // <li> the minimum nll, fitstatus and convergence quality for each fit </li> 
	     // <li> for each fit and for each non-constant parameter, the value, error and pull of the parameter are stored </li>
	     // </ul>
	     return fDetailedOutput;
     }
         
     // set the conditional observables which will be used when creating the NLL
     // so the pdf's will not be normalized on the conditional observables when computing the NLL 
     virtual void SetConditionalObservables(const RooArgSet& set) {fConditionalObs.removeAll(); fConditionalObs.add(set);}

     virtual void SetVarName(const char* name) { fVarName = name; }
     virtual const TString GetVarName() const {return fVarName;}

     virtual RooAbsPdf * GetPdf() const { return fPdf; }

      
      //      const bool PValueIsRightTail(void) { return false; } // overwrites default

      void FillTree(Bool_t flag=true) {fFillTree=flag;}
      TTree* GetTree(Bool_t getOwnership=false);
      void FlagData(const RooAbsData* data) {fFlagData=data;}
      void SkipToys(Int_t flag=0) {fSkipToys=flag;}
      void SetComparisonTestStat(TestStatistic* ts) {fComparisonTestStat=ts;}
      void UseMinos(const RooArgSet* set=0, Bool_t flag=true) {fMinos=flag; fMinosSet=set;}
      void PlotParameterProfiles(Int_t flag=0) {fParmProf=flag;} // 0:just data, >0: also for a number of toys
      TList& GetInfo() {return fInfo;}
      void SaveToyData(RooWorkspace* ws,const RooArgSet* glob=0) {fSaveWS=ws; fSaveGlob=glob; fUseSaved=false;}
      void UseSavedToyData(RooWorkspace* ws) {fSaveWS=ws; fSaveGlob=0; fUseSaved=true;}
      ULong64_t SetOptimize(ULong64_t opt) {ULong64_t o=fOptimize; fOptimize=opt; return o;}
      ULong64_t SetOptimize(const char* opt=0, int verbose=-1);
      void SetInitialValues(const RooArgSet* init, const RooArgSet* initCond=0) {fInitialValues=init;fInitialValuesCond=initCond;}
      void SetFixInitially(const RooArgSet* fix) {fFixInitially=fix;}
      void SetHypoTestCalculator(RooStats::HypoTestCalculatorGeneric* calc) { fCalc=calc; }
      const RooFitResult* FitResult()     const { return resultU; }
      const RooFitResult* FitResultCond() const { return resultC; }
      RooAbsReal*   GetNll()       const { return fNll;       }
      RooMinimizer* GetMinimizer() const { return fMinim; }

      static TString Join(const RooAbsCollection& c, const TString& sep=",", bool showname=false);
      static RooArgSet* GetCommonSnapshot (const RooAbsCollection& a, const RooAbsCollection& b);
      static bool saveHessian (TList& info, const TString& hesseName, int verbose=1);
      TMatrixDSym reducedCovarianceMatrix (const RooFitResult* res, const RooArgList& params);
      void TCollectionRemoveAll (TCollection& c, const char* name);

      static const ULong64_t default_optimize;

  private:

      double GetMinNLL        (int stage, int& status, int& tries, int& nfcn, RooFitResult*& result,
                               RooArgSet* attachedSet, const RooArgSet* pois, const RooArgSet* reset=0);
      double GetMinNLL1       (int stage, int& status, int& tries, int& nfcn, RooFitResult*& result,
                               RooArgSet* attachedSet, const RooArgSet* pois, const RooArgSet* reset=0);
      double UnconditionalFit (int stage, int& status, int& tries, int& nfcn,
                               RooArgSet* attachedSet, const RooArgSet* pois);
      void PrintFitStatus (int stage, const char* msg, int status, const RooMinimizer& minim, const RooFitResult* res,
                           Double_t initialNll, Double_t initialPOI, const RooRealVar* poi,
                           bool verbose=true);

   private:

      RooAbsPdf* fPdf;
      RooAbsReal* fNll; //! pointer to negative log-likelihood function
      const RooArgSet* fCachedBestFitParams;
      const RooAbsData* fLastData;
      //      Double_t fLastMLE;
      LimitType fLimitType;
      Bool_t fSigned;
      
      // this will store a snapshot of the unconditional nuisance
      // parameter fit.
      int fDetailedOutputEnabled;
      bool fDetailedOutputWithErrorsAndPulls;
      RooArgSet* fDetailedOutput; //!
      RooArgSet fConditionalObs;    // conditional observables 
      
      TString fVarName;

      static Bool_t fgAlwaysReuseNll ;
      Bool_t fReuseNll ;
      TString fMinimizer;
      Int_t fStrategy;
      Double_t fTolerance, fTolerance0;
      Int_t fPrintLevel;

      Bool_t fFillTree;
      TTree* fTree;
      struct TreeVars* fVars; //!
      const RooAbsData* fFlagData; //!
      TestStatistic* fComparisonTestStat;
      Bool_t fMinos;
      Int_t fParmProf;
      TList fInfo;
      Int_t fToyNum;
      const RooArgSet* fMinosSet;
      Int_t fSkipToys, fSeed;
      Bool_t fUseSaved;
      RooWorkspace* fSaveWS;
      const RooArgSet* fSaveGlob;
      Int_t fProblems;
      ULong64_t fOptimize;
      const RooArgSet *fInitialValues, *fInitialValuesCond;
      Bool_t fSavedInitialSnapshot;
      Int_t fEvalNum;
      RooAbsReal* fNll0; //!
      const RooArgSet* fFixInitially;
      RooStats::HypoTestCalculatorGeneric* fCalc;
      RooFitResult *resultU, *resultC;
      RooMinimizer* fMinim;
      RooArgSet* fLastInitalSnap[3]; //!
      Bool_t fUseNll; //!

#ifndef __ACLIC__
   protected:

      ClassDef(ProfileLikelihoodTestStatEnhanced,22)   // implements the enhanced profile likelihood ratio as a test statistic to be used with several tools
#endif
   };
}

#ifdef __ACLIC__
#include "ProfileLikelihoodTestStatEnhanced.cxx"
#endif

#endif
