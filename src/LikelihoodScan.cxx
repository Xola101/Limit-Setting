/*==============================================================================
$Id: LikelihoodScan.C 515454 2012-08-29 20:14:34Z adye $

Likelihood plots.

Usage from the command-line (if an executable has been build using "make"):

  LikelihoodScan WORKSPACE-FILE.root \
      -w workspace -m modelConfig -d dataset \
      -p jobNum -t nToys -M invMass -0 muMin -1 muMax \
      -W nProofWorkers -C -O OPT-LEVEL -H -r OUTPUT-RESULT-FILE

The calculation is performed by the external class, LikelihoodCalculator.

==============================================================================*/

#define ALLOW_INCLUDE_CODE
#include "LikelihoodCalculator.h"

#if defined(__CINT__) || defined(__ACLIC__)

#include "TROOT.h"

LikelihoodCalculator* calc= 0;

static void ResetCalc()
{
  if (!calc) return;
  // If run interactively, remove canvas and all histograms that might have been
  // created with a previous invocation.
  delete calc; calc= 0;
  gDirectory->Clear();
}

void LikelihoodScan (const char* args="")
{
  // Main routine when run from the ROOT prompt
  ResetCalc();
  calc= new LikelihoodCalculator ("LikelihoodScan", args);
  calc->Run();
}

#else

int main (int argc, const char** argv)
{
  // Main program when run stand-alone
  LikelihoodCalculator calc ("LikelihoodScan", argc, argv);
  return calc.Run();
}
#endif
