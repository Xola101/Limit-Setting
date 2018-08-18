/*==============================================================================
$Id: StandardFrequentistDiscovery.C 515454 2012-08-29 20:14:34Z adye $

Cross-section significance calculation using toys.
This version allows for the work to be split up into many batch or Grid jobs.
Uses RooStats::FrequentistCalculator to compute the p-value.

Usage from the command-line (if an executable has been build using "make"):

  StandardFrequentistDiscovery WORKSPACE-FILE.root \
      -w workspace -m modelConfig -d dataset \
      -p jobNum -t nToys -M invMass -0 muMin -1 muMax \
      -W nProofWorkers -C -O OPT-LEVEL -H -r OUTPUT-RESULT-FILE

or from the ROOT prompt (change other options in the parameter defaults at the top of the file):

  root [0] .x StandardFrequentistDiscovery.C("WORKSPACE-FILE.root","combWS","ModelConfig","combData")

The calculation is performed by the external class, PValueCalculator.

Author: Tim Adye, based on $ROOTSYS/tutorials/roostats/StandardFrequentistDiscovery.C from ROOT 5.34.00.

==============================================================================*/

#define ALLOW_INCLUDE_CODE
#include "PValueCalculator.h"

#if defined(__CINT__) || defined(__ACLIC__)

#include "TROOT.h"

PValueCalculator* calc= 0;

static void ResetCalc()
{
  if (!calc) return;
  // If run interactively, remove canvas and all histograms that might have been
  // created with a previous invocation.
  delete calc; calc= 0;
  gDirectory->Clear();
}

double StandardFrequentistDiscovery(
   const char* infile,
   const char* workspaceName,
   const char* modelConfigNameSB = "ModelConfig",
   const char* dataName = "combData",
   int toys = 1000,
   double poiValueForBackground = 0.0,
   double poiValueForSignal = 1.0
) {
  // StandardFrequentistDiscovery with function arguments as in $ROOTSYS/tutorials/roostats/StandardFrequentistDiscovery.C
  ResetCalc();
  calc= new PValueCalculator ("StandardFrequentistDiscovery");
  calc->fileName.push_back(infile);
  calc->wsName= workspaceName;
  calc->modelSBName= modelConfigNameSB;
  calc->dataName= dataName;
  calc->ntoys = toys;
  calc-> bPOI.assign (1, poiValueForBackground);
  calc->sbPOI.assign (1, poiValueForSignal);
  calc->Run();
  if (!calc->res) return -1.0;
  return calc->res->NullPValue();
}


void StandardFrequentistDiscovery (const char* args="")
{
  // Main routine when run from the ROOT prompt
  ResetCalc();
  calc= new PValueCalculator ("StandardFrequentistDiscovery", args);
  calc->Run();
}

#else

int main (int argc, const char** argv)
{
  // Main program when run stand-alone
  PValueCalculator calc ("StandardFrequentistDiscovery", argc, argv);
  return calc.Run();
}
#endif
