//=====================================================================-*-C++-*-
// File: $Id: LikelihoodCalculator.h 600726 2014-06-07 01:58:07Z adye $
//==============================================================================

#ifndef LikelihoodCalculator_h
#define LikelihoodCalculator_h

#include "WorkspaceCalculator.h"

#if defined(ALLOW_INCLUDE_CODE) && defined(__CINT__) && !defined(__ACLIC__)
static int compiled_LikelihoodCalculator = gSystem->CompileMacro("LikelihoodCalculator.cxx", "k");
#else

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
//#define USE_THn
#endif

#include "TVectorDfwd.h"
#include "TMatrixDfwd.h"
class TH1;
class TObjArray;
class TGraph;
class TAxis;

#ifdef USE_THn
class THnBase;
template <class T> class THnT;
typedef THnT<Double_t> THnD;
#else
class TArrayD;
class THnSparse;
template <class CONT> class THnSparseT;
typedef THnSparseT<TArrayD> THnSparseD;
#endif

class LikelihoodCalculator : public WorkspaceCalculator {
  public:
#ifdef USE_THn
  typedef THnBase    HistTypeBase;
  typedef THnD       HistType;
#else
  typedef THnSparse  HistTypeBase;
  typedef THnSparseD HistType;
#endif

//==============================================================================
// Constructors and destructor
//==============================================================================

  LikelihoodCalculator (const char* name="LikelihoodCalculator");
  LikelihoodCalculator (const char* name, int argc, const char* const* argv);
  LikelihoodCalculator (const char* name, const char* args);
  void Init();
  virtual ~LikelihoodCalculator();

//==============================================================================
// Control methods
//==============================================================================

  using WorkspaceCalculator::ParseArgs;
  virtual int ParseArgs (int argc, const char* const* argv);
  virtual void Usage() const;
  virtual int RunOverWorkspace();
  virtual int RunCalculator();
  virtual void AdjustResults (bool readResult);
  virtual void SetDefaults();
  virtual void ShowParms() const;

//==============================================================================
// Other methods
//==============================================================================

  // get limits, print them, and write them to an ntuple
  virtual TTree* GetLimits();
  // read the results from files
  using WorkspaceCalculator::ReadResult;
  virtual int ReadResult (TFile* f, const char* objname);
  // print a summary from the results
  virtual void PrintResult();
  // plot the results
  using WorkspaceCalculator::PlotResult;
  virtual void PlotResult(TPad* canvas);

  virtual int ReadTree (TFile* f, UInt_t& thisSeed);

  virtual bool BkgIsAlt() { return true; };

  void PlotProfile (TPad* canvas, const TH1D* h, const TMatrixD& xval, Double_t mle=NaN, int axes=0);
  Double_t Plot2D (TPad* canvas, const HistTypeBase* h, const std::vector<Int_t>& axes= std::vector<Int_t>(), TMatrixD* path=0);
  void DrawPlots (TPad* canvas, const TCollection& plots);
  void UseFittedPath (TMatrixD& path) const;
  TH1D* Profile2D (const HistTypeBase* h, int axis,
                   TMatrixD* path, const TMatrixD* oldPath, const std::vector<Int_t>& axes,
                   Double_t& mle, Double_t& xmin, Double_t& ymin, TMatrixD& xval);
  void FixedErrors (const TH1D* h, const TMatrixD& xval, const TObjArray& hist1d);

//==============================================================================
// Utility routines
//==============================================================================

  static HistTypeBase* ReadResult (TFile* f, const char* objname, int verbose);
  static int MergeResult (HistTypeBase*& res, const HistTypeBase* radd, const char* fname);
  static std::vector<int> GetBin (int ibin, const std::vector<int>& npoints, bool lr=true, int offset=0, int dfix=-1, int vfix=0);
  static int GetBin (const std::vector<int>& bin, const std::vector<int>& npoints, bool lr=true, int offset=0, int dfix=-1);
  static TString BinName (const HistTypeBase* h, const std::vector<Int_t>& bin, Int_t skip=-1);
  static TString HistShape (const HistTypeBase* h);
  static void InterpolateMinimum (double xl, double x0, double xh, double yl, double y0, double yh,
                                  double &xm, double &ym);
  static Int_t FindErrors (const TH1D* h, Double_t xmin, const std::vector<Double_t>& up, TMatrixD& xval,
                           int interpolate=0, int verbose=0);
  static void FindExtrema (const HistTypeBase* h, Double_t mle=0.0, int interpolate=1, int verbose=0, int nup=2);

  static HistTypeBase* Profile (const HistTypeBase* h, Int_t axis, TObjArray* h1d=0,
                                TMatrixD* path=0, const TMatrixD* oldPath=0, const std::vector<Int_t>& axes=std::vector<Int_t>(),
                                int interpolate=1, int verbose=0);
  static Double_t FindMinimum (const TH1D* h, Int_t& imin, Double_t& xmin, int interpolate=0, int verbose=0);
  static void HistToLikelihood (TH1* h, Double_t mle);
  static std::vector<Double_t> GetUpLevels (int nupLevels, int ndim=1, int verbose=0, bool lim95=false);
  static TGraph* Path2Graph (const TMatrixD& path, int axis, const std::vector<Int_t>& axes);
  static TAxis* ReplaceAxis (TAxis* dest, const TAxis* source);
  static std::vector<Double_t> GetBinPoints (const TAxis* ax);
  static Double_t GetBinPoint (const TAxis* ax, Int_t bin);
  static Double_t GetBinPoint (const HistTypeBase* h, Int_t d, Int_t bin);
  static Int_t FindBinExt (const TAxis* ax, Double_t x, bool upper=false);
  static Double_t GetBinLowEdgeExt (const TAxis* ax, Int_t bin);
  static void PrintErrorPoints (const TH1D* h, const TMatrixD& xval, TMatrixD* path, const std::vector<Double_t>& upLevel);

//==============================================================================
// Private utility routines
//==============================================================================
private:
  int AddMassAxis  (HistTypeBase*& res, std::vector<HistTypeBase*>& res2) const;
  int ResizeResult (HistTypeBase*& res, std::vector<HistTypeBase*>& res2) const;

//==============================================================================
// Data members (public for simple access)
//==============================================================================
public:
  std::vector<int> npoints;
  std::vector<bool> fitInRange;
  int         nscan, nscanDefault;
  int         npoi;
  std::vector<Double_t> scanMin;
  std::vector<Double_t> scanMax;
  int         firstPoint;
  int         lastPoint;       // Set >=0 to scan only a subset of the points
  std::vector< std::vector<Double_t> > pointsToScan;
  int         jobSet;
  int         interpolate;
  int         plotAllProjections;
  TList       toPlot;
  HistTypeBase* res;
  std::vector<HistTypeBase*> res2;
  std::vector<Double_t> markerSM;
  Double_t    ymax;
  std::vector<std::string> extraPOI;
  int         ncontour2D, nupLevels;
  bool        lim95;
  int         massAxis;
  TString     plotOpt;
  std::vector<Double_t> upLevel, upLevel2D;

#if !defined(ALLOW_INCLUDE_CODE) || !defined(__ACLIC__)
protected:
  ClassDef(LikelihoodCalculator,3)    // implements the likelihood scan
#endif
};

#if defined(ALLOW_INCLUDE_CODE) && defined(__ACLIC__)
#include "LikelihoodCalculator.cxx"
#endif

#endif
#endif
