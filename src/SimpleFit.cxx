/*==============================================================================
$Id: SimpleFit.C 550136 2013-06-06 19:01:21Z adye $

Simple fit to data in workspace.

Usage from the command-line (if an executable has been build using "make"):

  SimpleFit WORKSPACE-FILE.root \
      -w workspace -m modelConfig -d dataset \
      -W nProofWorkers -O OPT-LEVEL -H -r OUTPUT-RESULT-FILE

The calculation is performed by the external class, FitCalculator.

==============================================================================*/

#define ALLOW_INCLUDE_CODE
#include "FitCalculator.h"

#if defined(__CINT__) || defined(__ACLIC__)

#include "TROOT.h"

FitCalculator* calc= 0;

static void ResetCalc()
{
  if (!calc) return;
  // If run interactively, remove canvas and all histograms that might have been
  // created with a previous invocation.
  delete calc; calc= 0;
  gDirectory->Clear();
}

void SimpleFit (const char* args="")
{
  // Main routine when run from the ROOT prompt
  ResetCalc();
  calc= new FitCalculator ("SimpleFit", args);
  calc->Run();
}

#else
  printf("StandAlone mode disabled\n");
//int main (int argc, const char** argv)
//{
  // Main program when run stand-alone
//  FitCalculator calc ("SimpleFit", argc, argv);
//  return calc.Run();
//}
#endif
