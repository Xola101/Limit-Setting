//=====================================================================-*-C++-*-
// File: $Id: PValueCalculator.h 548359 2013-05-21 22:04:57Z adye $
//==============================================================================

#ifndef PValueCalculator_h
#define PValueCalculator_h

#include "WorkspaceCalculator.h"

#if defined(ALLOW_INCLUDE_CODE) && defined(__CINT__) && !defined(__ACLIC__)
static int compiled_PValueCalculator = gSystem->CompileMacro("PValueCalculator.cxx", "k");
#else

class TH1;
class TH1D;
class TH2D;
class TF1;

namespace RooStats {
  class HypoTestResult;
}

class PValueCalculator : public WorkspaceCalculator {
  public:

//==============================================================================
// Constructors and destructor
//==============================================================================

  PValueCalculator (const char* name="PValueCalculator");
  PValueCalculator (const char* name, int argc, const char* const* argv);
  PValueCalculator (const char* name, const char* args);
  void Init();
  virtual ~PValueCalculator();

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
  virtual void PrintResult() { if (res) PrintResult(res); }
  virtual void PrintResult (RooStats::HypoTestResult* r);
  // plot the results
  using WorkspaceCalculator::PlotResult;
  virtual void PlotResult(TPad* canvas);

  virtual int ReadTree (TFile* f, UInt_t& thisSeed);

  using WorkspaceCalculator::DropBadSamples;
  virtual int DropBadSamples();
  using WorkspaceCalculator::ApplyCuts;
  virtual int ApplyCuts();
  // limit number of toys in result
  using WorkspaceCalculator::LimitSamples;
  virtual int LimitSamples();

  virtual void PlotTSdiff (TPad* canvas, TH1D* he, int tsType=2, Double_t obsTS=NaN, const TF1* fit=0);
  virtual void PlotPValue (TPad* canvas, TH1* h, RooStats::HypoTestResult* r, int tsType=2, const TF1* fit=0);
  virtual int TopCorr (const std::vector<std::string>& npnames,
                       const int npoi,
                       std::vector<int>& itop,
                       std::vector<int>& jtop,
                       std::vector<Double_t>& topcor,
                       TH2D** hcor_ret=0,
                       const Double_t maxCorrPOI=-1.0,
                       const Double_t maxCorr=-1.0,
                       const int maxTop=-1);
  virtual void PlotTS2 (TPad* canvas, int tsType, Double_t tsdata, Double_t muhat,
                        const char* poiname, const char* hname, const char* hename);
  virtual void PlotPOIvsTS (TPad* canvas, Double_t tsdata, const char* poiname);
  virtual void PlotVar (TPad* canvas, const char* name, const char* dname, const char* varcut, Long64_t idata, Double_t genval=NaN);

  virtual bool BkgIsAlt() { return false; };

//==============================================================================
// Utility routines
//==============================================================================

  // estimate errors using resampling
  static Double_t ResampleErrors (const RooStats::HypoTestResult* r, Double_t* experr, int ntoys=50);
  using WorkspaceCalculator::ResampleToy;
  static RooStats::HypoTestResult* ResampleToy (const RooStats::HypoTestResult* rold);
  static Double_t ExpectedPValue (RooStats::HypoTestResult* r, Double_t* experr=0, Double_t* q0exp=0, double nsig=0.0);
  static Double_t GetVal (TTree* tree, const char* name, Long64_t entry=-1, Double_t defval=NaN);

//==============================================================================
// Private utility routines
//==============================================================================
protected:

//==============================================================================
// Data members (public for simple access)
//==============================================================================
public:
  bool fitLeadbetter;             // fit Leadbetter term to q0 distribution
  RooStats::HypoTestResult* res;  // toy results

#if !defined(ALLOW_INCLUDE_CODE) || !defined(__ACLIC__)
protected:
  ClassDef(PValueCalculator,1)    // implements the p-value calculator
#endif
};

#if defined(ALLOW_INCLUDE_CODE) && defined(__ACLIC__)
#include "PValueCalculator.cxx"
#endif

#endif
#endif
