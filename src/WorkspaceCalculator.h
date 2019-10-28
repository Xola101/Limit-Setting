//=====================================================================-*-C++-*-
// File: $Id: WorkspaceCalculator.h 727292 2016-03-02 17:22:07Z adye $
//==============================================================================

#ifndef WorkspaceCalculator_h
#define WorkspaceCalculator_h

#if defined(ALLOW_INCLUDE_CODE) && defined(__CINT__) && !defined(__ACLIC__)
static int compiled_WorkspaceCalculator = gSystem->CompileMacro("WorkspaceCalculator.cxx", "k");
#else

#include <string>
#include <cstdlib>
#include <vector>
#include "Rtypes.h"
#include "TNamed.h"
#include "TString.h"
#include "TVectorDfwd.h"
#include "TList.h"
#include "TTime.h"
#include "TClass.h"
#include "TTree.h"
#include "TBranchElement.h"
#include "RooArgSet.h"

class TFile;
class TH1;
class TH1D;
class TPad;
class TF1;
class TGraph;
class TPaveText;

class RooWorkspace;
class RooArgList;
class RooRealVar;
class RooAbsCollection;
class RooDataSet;
class RooAbsData;
class RooAbsPdf;

namespace RooStats {
  class ModelConfig;
  class HypoTestResult;
  class SamplingDistribution;
  class TestStatistic;
  class ToyMCSampler;
  class HypoTestCalculatorGeneric;
  class ProfileLikelihoodTestStatEnhanced;
};

#ifndef DEFINED_NAN
#define DEFINED_NAN
#ifdef __CINT__
#include "TMath.h"
static const Double_t NaN = TMath::QuietNaN();
#else
#include <limits>
static const Double_t NaN = std::numeric_limits<Double_t>::quiet_NaN();
#endif
#endif

struct TOwnedList : public TList {
  // A collection class for keeping TObjects for deletion.
  // TOwnedList is like TList with SetOwner(), but really deletes all objects, whether or not on heap.
  // This is a horrible hack to work round the fact that RooArgSet and RooDataSet objects have have IsOnHeap() false.
  TOwnedList();
  virtual ~TOwnedList();
  virtual void Clear (Option_t* option="");
#if !defined(ALLOW_INCLUDE_CODE) || !defined(__ACLIC__)
  ClassDef(TOwnedList,0)
#endif
};

class WorkspaceCalculator : public TNamed {
  public:

  enum OptimizeBits {
    kReuseNLL            = 1LL <<  0,     // R: reuse NLL in test statistic
    kUseMultiGen         = 1LL <<  1,     // M: use multigen optimisation of toy MC generation
    kMinuit2             = 1LL <<  2,     // F: use Minuit2 (default is Minuit)
    kStrategy0           = 1LL <<  3,     // S: use minimizer strategy 0
    kVectorStore         = 1LL <<  4,     // V: use vector store in datasets
    kGenerateBinned      = 1LL <<  5,     // B: use binned toy generation for gg channel
    kInitialFit          = 1LL <<  6,     // I: perform initial fit
    kAdjustRanges        = 1LL <<  7,     // A: adjust parameter ranges
    kSetConst            = 1LL <<  8,     // C: set nuisance parameters with zero-range to constant
    kSkipDataTS          = 1LL <<  9,     // D: skip TS calculation for data
    kUseCondSnapshot     = 1LL << 10,     // s: use conditionalNuis_0/1 snapshots to generate toys instead of profiling
    kLoadInitialSnapshot = 1LL << 11,     // i: load initial snapshot 'ucmles' (initialises data fit and profiling)
    kMinosData           = 1LL << 12,     // m: run MINOS on conditional data fit
    kInterpolateErrors   = 1LL << 13,     // e: set error interpolation code 4 for all functions
    kSaveInfo            = 1LL << 14,     // r: save extra info in results file: random number seeds and fit results for each toy
    kUseSavedToys        = 1LL << 15,     // u: use toys saved in the workspace
    kStrategy2           = 1LL << 16,     // 2: use minimizer strategy 2
    kChainNtuple         = 1LL << 17,     // c: use chain instead of merging statsTree files (less mem/faster, but very slow with -H)
    kAddSLRTS            = 1LL << 18,     // L: run SimpleLikelihoodRatioTestStat as a second test stat
    kAnalytic            = 1LL << 19,     // a: analytic statistical errors in PDF
    kPLTSEnhanced        = 1LL << 20,     // E: Use ProfileLikelihoodTestStatEnhanced instead of RooStats::ProfileLikelihoodTestStat
    kInitGlobObs         = 1LL << 21,     // g: Use global observables to initialise fit
    kExpandRanges        = 1LL << 22,     // G: expand gg slope parameter range
    kDisableNumInt       = 1LL << 23,     // N: disable numeric integration for llll PDFs
    kFixStatError        = 1LL << 24,     // f: fix small histogram bin statistics nuisance parameters
    kFixInitially        = 1LL << 25,     // x: perform initial fit with constrained NPs fixed
    kReuseAltToys        = 1LL << 26,     // t: re-use the same toys for the alternate hypothesis
    kSaveWorkspace       = 1LL << 27,     // W: save workspace with fitted snapshots
    kNoFits              = 1LL << 28,     // X: skip all fits, just check setup
    kAddPseudoData       = 1LL << 29,     // P: add pseudo-data to 4l datasets
    kBinnedLikelihoodOpt = 1LL << 30,     // b: activate binned likelihood calculation for binned models
    kSetGlobalObservable = 1LL << 31,     // o: identify global observables for RooFit with GLOBAL_OBSERVABLE tag
    kFixCache            = 1LL << 32,     // p: fix cache for RooStarMomentMorph pdfs
    kAltMinos            = 1LL << 33,     // v: use alternate method for MINOS fits
    kFastInit            = 1LL << 34,     // O: skip slower init steps: don't print NPs split into constrained, unconstrainted, and gamma
    kAddErrors           = 1LL << 35,     // d: add nominal errors to NPs if not already set
    kUpdGlobalObservables= 1LL << 36,     // U: when adding or removing NPs, also add or remove their global observables
  };

//==============================================================================
// Constructors and destructor
//==============================================================================

  WorkspaceCalculator (const char* name="WorkspaceCalculator");
  void Init();
  void InitCalculator();
  virtual ~WorkspaceCalculator();
  virtual void DeleteCalculator();

//==============================================================================
// Control methods
//==============================================================================

  virtual int Run();
  virtual int RunOverWorkspace() = 0;
  virtual int WriteResultFile();
  virtual int SetupWorkspace();
  virtual int SetupMinimizer();
  virtual int SetupModel();
  virtual int SetupInitialParms();
  virtual int SetupTestStat();
  virtual int SetupSampler();
  virtual int SetupCalculator();
  virtual int RunCalculatorProof();
  virtual int RunCalculator() = 0;
  virtual int GetTestStatInfo();
  virtual void DoInitialFit();
  virtual void AdjustResults (bool readResult);
  virtual void AnalyseResults();

  virtual void ShowParms() const {};
  virtual int ParseArgs (const char* args="");
  virtual int ParseArgs (int argc, const char* const* argv) = 0;

//==============================================================================
// Other methods
//==============================================================================

  // get limits, print them, and write them to an ntuple
  virtual TTree* GetLimits() = 0;
  // print summary of SamplingDistribution
  virtual void PrintSamplingDistribution (RooStats::SamplingDistribution* s, const char* name, int nprint= 10);
  // drop bad samples (enable with dropNegative or dropNaN)
  virtual int DropBadSamples() {return 0;};
  virtual RooStats::SamplingDistribution* DropBadSamples (RooStats::SamplingDistribution* s);
  virtual int ApplyCuts() {return 0;};
  virtual RooStats::SamplingDistribution* ApplyCuts (RooStats::SamplingDistribution* s, bool use_alt=false);
  virtual int NumFailures (const RooStats::SamplingDistribution* s);

  virtual void SetDefaults();

  void ModifyInterpolationForAll(RooWorkspace* ws, int codeNorm, int codeShape);

  // read the results from files
  virtual int ReadResult (const std::vector<std::string>& fileNames);
  virtual int ReadResult (TFile* f, const char* objname) = 0;
  // print a summary from the results
  virtual void PrintResult() = 0;
  // limit number of toys in result
  virtual int LimitSamples() {return 0;};
  // plot the results
  virtual void PlotResult();
  virtual void PlotResult(TPad* canvas) = 0;

  virtual int ReadTree (TFile* f, UInt_t& thisSeed) = 0;

  virtual void wstitle();
  virtual void PrintPlot (TPad* canvas);
  virtual TString tsName (int tsType=-1);
  virtual TString chi2Name (int tsType=2, int ndf=-1);
  virtual void PlotTS (TPad* canvas, TH1D* h, TH1D* ha, TH1D* hb, int tsType, Double_t obsTS, TF1** fitOut=0, bool norm=true, Color_t lineCol=kBlue+2);
  virtual TF1* AsymTF1 (int tsType, int typ=0, const TH1* h=0, Double_t norm=1.0);
  virtual bool BkgIsAlt() = 0;
  virtual int SetOptimize(const char* a, const char* opt="");
  virtual void PlotDataSet (TPad* canvas, const RooAbsPdf& pdf, const RooAbsData& data, bool drawComponent=true);

  virtual Double_t CondFit   (const RooRealVar* par, Double_t val, RooAbsCollection* fitted=0, const char* opt=0);
  virtual Double_t DoFit     (const RooArgSet& poiSet, int type=0, RooAbsCollection* fitted=0, const char* opt=0);
  virtual int      RunMinos  (const RooArgSet& pars, Double_t nll_min=NaN);
  virtual int      findSigma (Double_t nll_min, RooRealVar* par, Double_t nsigma=1.0,
                              RooAbsCollection* fitted_lo=0, RooAbsCollection* fitted_hi=0,
                              const char* opt=0, Double_t precision=-1.0);
  virtual  Double_t findSigma (Double_t nll_min, Double_t val_guess, Double_t val_mle, const RooRealVar* par,
                               Double_t nsigma=+1.0, RooAbsCollection* fitted=0, const char* opt=0, Double_t precision=-1.0);

  virtual RooArgSet* SetParameters (const RooAbsCollection* setvars, bool updateNP=true);

//==============================================================================
// Utility routines
//==============================================================================
  static void PrintOptimize (ULong64_t opt, const char* first="");
  static const RooArgSet* SetValueInSnapshot (RooStats::ModelConfig* mc, const std::vector<double>& val);

  static Bool_t LoadSnapshot (RooWorkspace* w, const TString& snapName, const char* msgfmt=0);
  static RooArgSet* GetSnapshot (RooWorkspace* w, const TString& snapName, const char* msgfmt,
                                 const RooArgSet* vars, const char* resetSnap="initialNuisAndPoi");

  static void AdjustParameterRanges (RooWorkspace* ws, const RooArgSet* parms,
                                     const RooAbsPdf* pdf, const RooArgSet* globalObservables,
                                     double stdDevRange=10.0, double widerRange=30.0,
                                     const char* snapshotName="conditionalNuis_muhat",
                                     int verbose=0, int adjustRanges=1, double fixStatError=-1.0);
  static int NominalValue (const RooWorkspace* ws, const TString& name, Double_t& nomVal, Double_t& nomErr);
  static RooArgSet* SelectConstrained (const RooArgSet* parms,
                                       const RooAbsPdf* pdf, const RooArgSet* globalObservables,
                                       bool constrained=true, int verbose=0);


  template <class T>
  static Int_t SetBranchAddress (TTree* nt, const char *bname, T *add, bool optional=false);
  template <class T>
  static Int_t SetBranchAddress (TTree* nt, const char *bname, std::vector<T>* vec, bool optional=false, const T& init=T());

  // limit number of toys in result
  static RooStats::SamplingDistribution* LimitSamples (RooStats::SamplingDistribution* s, int ntoys, int skip=0);
  // estimate errors using resampling
  static RooStats::SamplingDistribution* ResampleToy (RooStats::SamplingDistribution* s);

  static Double_t FillSamples (TH1* h, const RooStats::SamplingDistribution* s, bool shiftZero= false);
  static TH1* NewShaded (const TH1* h, Double_t minShaded, Style_t fillStyle=3004, Color_t fillColor=-1);
  static TH1D* NewTSHist (const char* name, const char* title, double ref=NaN, bool twoSided=false, int plotType=0,
                          double negVal=-0.025);
  static Long64_t Project (TTree* nt, TH1* h, const char* varexp, const char* selection="",
                           Option_t* option="", Long64_t nentries=1000000000, Long64_t firstentry=0);

  static void PrintResourcesUsed (const TTime& start, const char* msg=0);
  static void PrintResourcesUsed (const char* msg=0);
  static Int_t ReadObject (const char* name, TObject* obj, TDirectory* dir=0);

  static TString Join (const std::vector<int>&    v, const TString& sep=",");
  static TString Join (const std::vector<double>& v, const TString& sep=",", bool trim=true);
  static TString Join (const TVectorD&            v, const TString& sep=",", bool trim=true);
  static TString Join (const std::vector<double>& v1, const std::vector<double>& v2,
                       const TString& sep1=":", const TString& sep2=",", bool trim=true);
  static TString Join (const std::vector<std::string>& v, const TString& sep=",", bool trim=true);
  static TString Join (const RooAbsCollection& c, const TString& sep=",", int show=0,
                       const TString& sep0="=", const TString& sep1=":", const TString& flagc="C", const TString& flagu="-");

  static const char* Scan (const char* a, std::vector<int>& v, char sep=',', int nullVal=0, bool append=false);
  static const char* Scan (const char* a, std::vector<double>& v1, std::vector<double>& v2,
                           char sep1=':', char sep2=',', bool ok1=false, bool append=false);
  static const char* Scan (const char* a, RooAbsCollection& vars,
                           char sep0='=', char sep1=':', char sep2=',', char flagc='C', char flagn='-', char inf='-', bool append=true);
  static const char* Scan (const char* a, std::vector<double>&      v, char sep=',', bool append=false);
  static const char* Scan (const char* a, std::vector<std::string>& v, char sep=',', bool append=false);

  static size_t ReadFile (const char* filename, std::vector<std::string>& v, bool strip=true, char comment='#', bool append=true);
  static std::string ReadFile (const char* filename);

  static RooDataSet* SetDatasetBinning (const RooAbsPdf* pdf, const RooAbsData* data,
                                        Int_t setNbins=500, const char* generateBinnedTag=0,
                                        const char*   binnedCategories=0,
                                        const char* unbinnedCategories=0,
                                        int verbose=0,
                                        const char* weightVarName= "weightVar");
  static int AddPseudoData (RooAbsData* data, const RooAbsPdf* pdf, const char* categories, const char* binning,
                            int verbose=0, Double_t eps=1e-10);
  static Double_t AsymFunc     (const Double_t* x, const Double_t* p);
  static Double_t AsymPValue   (const Double_t* x, const Double_t* p);
  static Double_t AsymQuantile (const Double_t* x, const Double_t* p);
  static TGraph*  PValueGraph  (TF1* f);
  static TPaveText* AtlasLabel (TPad* canvas, Double_t xoff=0.05, Double_t yoff=0.0, bool transparent=false);
  static Int_t IncludeOverflows (TH1* h);
  static Int_t RemoveConstantParameters (RooAbsCollection* coll, const RooWorkspace* ws=0, const char* msg=0);
  static Double_t UseLimits (const RooRealVar* par, Double_t val, int verbose=1);
  static RooArgSet* getAllConstraints (const RooAbsPdf& pdf, const RooArgSet& observables, const RooArgSet* constrainedParams=0, int verbose=0);
  static RooArgSet* findGlobalObservable (const RooAbsPdf* pdf, const RooArgSet* observables, const RooArgSet* np=0, const RooArgSet* globalObservables=0,
                                          const RooArgSet* constraints=0, int verbose=0);
  static RooAbsArg* client (const RooAbsArg* a, const TClass* cls1=0, const TClass* cls2=0);
  static RooArgSet servers (const RooAbsArg* a);
  static bool CheckGaussian (RooRealVar* v, const RooArgSet* gobs, bool verbose=true);
  static bool CheckPoisson (RooRealVar* v, bool verbose=true);
  static bool CheckGamma (RooRealVar* v, bool verbose=true);
  static void AddErrors (RooStats::ModelConfig* mc, bool verbose=true);

  

//==============================================================================
// Private utility routines
//==============================================================================
protected:
  template <class PLTS>
  PLTS* SetupPLTS (PLTS* profll, const char* tsname= "profile likelihood test statistic");
  void SetInitialSnapshotTS (RooStats::ProfileLikelihoodTestStatEnhanced* profll);
  void SetInitGlobObs (RooStats::ProfileLikelihoodTestStatEnhanced* profll);
  void AddWeights (TTree*& tree, RooStats::HypoTestResult* r);
  Int_t GetWeights (RooAbsData* ds, const RooStats::SamplingDistribution* samp, Int_t is_alt,
                    std::vector<Double_t>& tsval, std::vector<Double_t>& wgt,
                    std::vector<Int_t>& alt, std::vector<Int_t>& label);

//==============================================================================
// Data members (public for simple access)
//==============================================================================
public:

  static const ULong64_t default_optimize;

  ULong64_t   optimize;
  int         plotResult;            // plot results? 0=no, 1=with captions, 2=without captions
  bool        writeResult;           // write result to a file
  int         nworkers;              // set non-zero to use PROOF Light when using toys (for freq or hybrid)
  bool        noSystematics;         // force all systematics to be off (i.e. set all nuisance parameters as constat
                                     // to their nominal values)
  int         initialFit;            // do a first  fit to the model (-1 : default, 0 skip fit, 1 always do fit)
  int         ntoysAlt;              // number of alt toys (defaults to 0 or 1)
  double      nToysRatio;            // ratio Ntoys S+B/B. Default is for no alt toys.
  int         nToysResample;         // number of toys to use in resampling error calculation
  bool        nToysLimit;            // limit number of toys in read sample
  bool        dropNegative;          // drop negative test statistic values?
  bool        dropNaN;               // drop NaN test statistic values (can occur with ROOT 5.33.02 ToyMCSampler)?
  double      negVal;                // actually only count values less than this as negative
  double      dropLarge;             // Drop TS values larger than this
  double      invMass;               // Only required for final result ntuple
  double      fitTol;                // Migrad tolerance: target EDM=0.001*fitTol.
                                     // RooFit default is 1, giving convergence when EDM<0.001.
  double      fitTol0;               // Initial fit with larger tolerance.
  double      fitPrec;               // Machine precision for Minuit, relative to Minuit2 default (~8.9e-16)
  double      OneSidedPositiveMin;   // Use this poimin with ProfileLikelihoodTestStat::SetOneSidedPositive (testStatType=5)
  TString     minimizerType;         // minimizer type (default is Minuit2 (kMinuit2) or Minuit)
  std::string resultFileName;
  TString     cut;
  TString     plotFile;
  int         verbose;
  int         doStatsTree;
  int         detailedOutput;        // save detailed output from test statistic (2=withErrorsAndPulls)
  int         skipPlot;
  bool        plotMuHatPositive;     // Plot NP only for muhat>=0
  std::vector<double>  bPOI;         // POI values to use for background
  double      bPOIdefault;
  std::vector<double> sbPOI;         // POI value to use for signal
  std::vector<double> nullPOI, altPOI;
  double      sbPOIdefault;
  bool        force;                 // don't exit on consistency error
  TString     initialSnapshot;       // loaded if optimize|kLoadInitialSnapshot (can be comma-sepaarted list)
  TString     initialSnapshotTS;     // initial snapshot before performing fits in ProfileLikelihoodTestStatEnhanced
  TString     nullMLESnapshot;       // use this snapshot instead of profiling null model
  TString     altMLESnapshot;        // use this snapshot instead of profiling alt model
  double      newParameterRanges;    // #sigma to use for parameter range adjustment
  int         parmProf;              // plot parameter profiles: 0=just data, >0: also for a number of toys
  int         maxFunctionCalls;      // RooMinimizer uses 500*N,
                                     // raw ROOT::Fit::Fitter/Minuit1/2 uses 200+100*N+5*N*N if 0 (N=num fit variables)
                                     // NB. RooMinimizer always uses 500*N regardless of what we set here - unless hacked
  int         minosSetNum;           // set number of parameters to determine using MINOS
  int         minosSetSize;          // number of paramters to run through MINOS at a time (0=all)
  std::vector<std::string> minosSetNames; // names  of paramters to run through MINOS
  RooArgSet*  minosSet;              // set of parameters to determine using MINOS
  int         skipToys;              // skip fit for this number of toys
  long int    seed;                  // use seed based on job num, else a unique (UUID) seed each job (see TRandom3::SetSeed)
  double      fixStatErrorMin;       // fix histogram bin statistics nuisance parameters when error is less than this

  std::string wsfile;
  TOwnedList  kept, saved, calcObjs, plotObjs; //! objects we need to keep. Note they are deleted in order: calcObjs, saved, kept.
  TTree*      statsTree;
  double      poihat, poierr;
  std::vector<double> poimin, poimax;
  std::vector<std::string> poiName;
  bool        doneRead;
  int         testStatType;          // test statistic type
  int         overrideTestStatType;  // >=0: multiple test stat types in different data files or specified by user
  int         samplerType;           // -1:no importance sampling, 0:adaptive importance sampling,
                                     // >0: importance sampling with this number of alt sets
  TString     optimizeTS;            // option flags to pass to ProfileLikelihoodTestStatEnhanced
  bool        showTitles, showCaptions;
  int         atlasStyle;
  int         haveWeights;
  double      nStdDevOverlap;        // overlap of importance samples
  RooWorkspace* wsInfo;

  std::vector<std::string> fileName;
  TString wsName;
  TString resultName;
  TString modelSBName;
  TString modelBName;
  TString dataName;
  int calculatorType;
  int ntoys;                         // number of signal or background toys
  int useNumberCounting;             // using number counting events? 0:no, 1:yes, -1:guess (default)
  TString nuisPriorName;             // name of prior for the nuisance. This is often expressed as constraint term in the global model.
                                     // It is needed only when using the HybridCalculator (calculatorType=1).
                                     // If not given, the prior PDF from ModelConfig is used.
  int initErr;
  TString seedName;
  int jobNum;
  RooWorkspace* ws;
  RooAbsData* data;
  RooStats::ModelConfig *sbModel, *bModel, *nullModel, *altModel;
  RooAbsPdf *sbPdf, *bPdf, *nullPdf, *altPdf;
  const RooArgSet *nullSnapshot, *altSnapshot;
  RooStats::TestStatistic *testStat, *testStat2;
  RooStats::ToyMCSampler* toymcs;
  RooStats::HypoTestCalculatorGeneric* calc;
  const char *nullModelName, *altModelName;
  std::string cmdArgs;
  RooArgList* varSettings;
  double asimovSig;
  int interpCodeNorm, interpCodeShape;
  std::vector<std::string> wsEdit;
  bool altResultName, showAtlasLabel, runLimit, runToys, runScan, constrainedOnly, optionalResult;
  RooArgSet *allParams, *nullParams, *minosParams;
  RooAbsCollection* minos_results;
  const char* resetGlobsSnapshot;

#if !defined(ALLOW_INCLUDE_CODE) || !defined(__ACLIC__)
protected:
  ClassDef(WorkspaceCalculator,11)    // implements the workspace calculator
#endif
};

template <class T>
Int_t WorkspaceCalculator::SetBranchAddress (TTree* nt, const char* bname, T* add, bool optional)
{
  // Shortcut for nt->SetBranchStatus(bname,1) and nt->SetBranchAddress(bname,add).
  UInt_t found=1;
  nt->SetBranchStatus (bname, 1, optional ? &found : 0);
  if (!found) return TTree::kMissingBranch;
  return nt->SetBranchAddress(bname,add);
}

template <class T>
Int_t WorkspaceCalculator::SetBranchAddress (TTree* nt, const char *bname, std::vector<T>* vec, bool optional, const T& init)
{
  // Shortcut for nt->SetBranchStatus(bname,1) and nt->SetBranchAddress(bname,vec).
  // If the branch is actually a scalar, not a vector, then just fill 1st element.
  UInt_t found=1;
  nt->SetBranchStatus (bname, 1, optional ? &found : 0);
  if (!found) return TTree::kMissingBranch;
  TBranchElement* be;
  if ((be= dynamic_cast<TBranchElement*>(nt->GetBranch(bname))) &&
      (be->GetClass() == TClass::GetClass(typeid(std::vector<T>)))) {
    vec->clear();
    be->SetObject(vec);
    return TTree::kMatch;
  } else {
    vec->resize(1,init);
    return nt->SetBranchAddress(bname,&vec->front());
  }
}

// Useful shortcuts for option parsing
inline long int Strtol (const char*& p, int base=10) { return strtol (p, const_cast<char**>(&p), base); }
inline double   Strtod (const char*& p)              { return strtod (p, const_cast<char**>(&p));       }

#if defined(ALLOW_INCLUDE_CODE) && defined(__ACLIC__)
#include "WorkspaceCalculator.cxx"
#endif

#endif
#endif
