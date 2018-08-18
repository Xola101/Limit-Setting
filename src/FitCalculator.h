//=====================================================================-*-C++-*-
// File: $Id: FitCalculator.h 582880 2014-02-12 20:06:15Z adye $
//==============================================================================

#ifndef FitCalculator_h
#define FitCalculator_h

#include "WorkspaceCalculator.h"

#if defined(ALLOW_INCLUDE_CODE) && defined(__CINT__) && !defined(__ACLIC__)
static int compiled_FitCalculator = gSystem->CompileMacro("FitCalculator.cxx", "k");
#else

class FitCalculator : public WorkspaceCalculator {
  public:

//==============================================================================
// Constructors and destructor
//==============================================================================

  FitCalculator (const char* name="FitCalculator");
  FitCalculator (const char* name, int argc, const char* const* argv);
  FitCalculator (const char* name, const char* args);
  void Init();
  virtual ~FitCalculator();

//==============================================================================
// Control methods
//==============================================================================

  using WorkspaceCalculator::ParseArgs;
  virtual int ParseArgs (int argc, const char* const* argv);
  virtual void Usage() const;
  virtual int RunOverWorkspace();
  virtual int RunCalculator();
//  virtual void AdjustResults (bool readResult);
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

  virtual bool BkgIsAlt() { return false; };

//==============================================================================
// Utility routines
//==============================================================================

//==============================================================================
// Private utility routines
//==============================================================================
private:

//==============================================================================
// Data members (public for simple access)
//==============================================================================
public:
  Double_t    nll;
  int         errorAnalysis;
  int         nParmsNow;
  Double_t    nsigma;
  bool        postfit;

#if !defined(ALLOW_INCLUDE_CODE) || !defined(__ACLIC__)
protected:
  ClassDef(FitCalculator,2)    // implements the likelihood scan
#endif
};

#if defined(ALLOW_INCLUDE_CODE) && defined(__ACLIC__)
#include "FitCalculator.cxx"
#endif

#endif
#endif
