//=====================================================================-*-C++-*-
// File: $Id: LimitCalculator.h 523421 2012-10-26 17:11:52Z adye $
//==============================================================================

#ifndef LimitCalculator_h
#define LimitCalculator_h

#include "WorkspaceCalculator.h"

#if defined(ALLOW_INCLUDE_CODE) && defined(__CINT__) && !defined(__ACLIC__)
static int compiled_LimitCalculator = gSystem->CompileMacro("LimitCalculator.cxx", "k");
#else

namespace RooStats {
  class HypoTestInverter;
  class HypoTestInverterResult;
}

class LimitCalculator : public WorkspaceCalculator {
  public:

//==============================================================================
// Constructors and destructor
//==============================================================================

  LimitCalculator (const char* name="LimitCalculator");
  LimitCalculator (const char* name, int argc, const char* const* argv);
  LimitCalculator (const char* name, const char* args);
  void InitCalculator();
  void Init();
  virtual ~LimitCalculator();
  virtual void DeleteCalculator();

//==============================================================================
// Control methods
//==============================================================================

  using WorkspaceCalculator::ParseArgs;
  virtual int ParseArgs (int argc, const char* const* argv);
  virtual void Usage() const;
  virtual int RunOverWorkspace();
  virtual int SetupScan();
  virtual int RunCalculator();
  virtual void AdjustResults (bool readResult);
  virtual int RebuildExpected();
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
  virtual void PrintResult() { PrintResult(res); }
  virtual void PrintResult (RooStats::HypoTestInverterResult* r);
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
  virtual bool BkgIsAlt() { return true; };

//==============================================================================
// Utility routines
//==============================================================================

  // estimate errors using resampling
  void ResampleErrors (const RooStats::HypoTestInverterResult* r, Double_t* obserr=0, Double_t* expserr=0, Double_t* expberr=0, int ntoys=50);
  using WorkspaceCalculator::ResampleToy;
  RooStats::HypoTestInverterResult* ResampleToy (const RooStats::HypoTestInverterResult* r);
  static Double_t ExpectedPValue (RooStats::HypoTestResult* r, Double_t* experr=0, Double_t* q0exp=0, double nsig=0.0);

//==============================================================================
// Private utility routines
//==============================================================================
protected:

//==============================================================================
// Data members (public for simple access)
//==============================================================================
public:
  RooStats::HypoTestInverterResult* res;
  bool        useCLs;
  int         npoints;
  double      scanMin;
  double      scanMax;
  double      confidenceLevel;
  int         rebuild;         // re-do extra toys for computing expected limits and rebuild test stat distributions
                               // N.B this requires much more CPU (factor is equivalent to nToyToRebuild)
  int         nToyToRebuild;   // number of toys used to rebuild
  int         firstPoint;
  int         lastPoint;       // Set >=0 to scan only a subset of the points
  bool        expCLsFromSB;    // expected CLs limits written to ntuple from CLs+b, else from background
  std::vector<Double_t> pointsToScan;
  int         jobSet;
  RooStats::HypoTestInverter* limitCalc;

#if !defined(ALLOW_INCLUDE_CODE) || !defined(__ACLIC__)
protected:
  ClassDef(LimitCalculator,1)    // implements the limit calculator
#endif
};

#if defined(ALLOW_INCLUDE_CODE) && defined(__ACLIC__)
#include "LimitCalculator.cxx"
#endif

#endif
#endif
