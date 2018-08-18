/*==============================================================================
$Id: StandardHypoTestInv.C 515454 2012-08-29 20:14:34Z adye $

Cross-section upper limit calculation using toys.
This version allows for the work to be split up into many batch or Grid jobs.
Uses RooStats::HypoTestInverter (inverted hypothesis test)
to perform a scan of the p-values and compute the limits.

Usage from the command-line (if an executable has been build using "make"):

  StandardHypoTestInv WORKSPACE-FILE.root \
      -w workspace -m modelConfig -d dataset \
      -a -p jobNum -n nPointsToScan -N points/job -P pointsToScan \
      -t nToys[,nToysBkg] -M invMass -0 scanMin -1 scanMax -2 poimax -c confidenceLevel \
      -W nProofWorkers -C -O OPT-LEVEL -H -r OUTPUT-RESULT-FILE

or from the ROOT prompt (change other options in the parameter defaults at the top of the file):

  root [0] .x StandardHypoTestInv.C("WORKSPACE-FILE.root","combWS","ModelConfig","combData")

The calculation is performed by the external class, LimitCalculator.

Author: Tim Adye, based on $ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C
                  from ROOT trunk revision 41812 (5.32.00-rc1+).

==============================================================================*/

#define ALLOW_INCLUDE_CODE

/* Project include */
#include "LimitCalculator.h"

#if defined(__CINT__) || defined(__ACLIC__)

/* Root framework include */
#include "TROOT.h"

LimitCalculator* calc= 0;

static void ResetCalc()
{
  if (!calc) return;
  // If run interactively, remove canvas and all histograms that might have been
  // created with a previous invocation.
  delete calc; calc= 0;
  gDirectory->Clear();
}

void
StandardHypoTestInv(const char * infile,
                    const char * wsName,
                    const char * modelSBName = "ModelConfig",
                    const char * modelBName = "",
                    const char * dataName = "combData",
                    int calculatorType = 0,
                    int testStatType = -1,  // default=2 or what's read from results file
                    bool useCLs = true ,
                    int npoints = 6,
                    double poimin =    0.0, // default=use workspace limits
                    double poimax = -999.0, // default=use workspace limits
                    int ntoys=1000,
                    bool useNumberCounting = false,
                    const char * nuisPriorName = 0){
  // StandardHypoTestInv with function arguments as in $ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C
  ResetCalc();
  calc= new LimitCalculator ("StandardHypoTestInv");
  calc->fileName.push_back(infile);
  calc->wsName= wsName;
  calc->modelSBName= modelSBName;
  calc->modelBName= modelBName;
  calc->dataName= dataName;
  calc->calculatorType = calculatorType;
  calc->testStatType = testStatType;
  calc->useCLs = useCLs;
  calc->npoints = npoints;
  if (poimin<poimax) {
    calc->poimin.assign (1, poimin);
    calc->poimax.assign (1, poimax);
  }
  calc->ntoys = ntoys;
  calc->useNumberCounting = useNumberCounting;
  if (nuisPriorName) calc->nuisPriorName = nuisPriorName;
  calc->Run();
}


void StandardHypoTestInv (const char* args="")
{
  // Main routine when run from the ROOT prompt
  ResetCalc();
  calc= new LimitCalculator ("StandardHypoTestInv", args);
  calc->Run();
}

#else

int main (int argc, const char** argv)
{
  // Main program when run stand-alone
  LimitCalculator calc ("StandardHypoTestInv", argc, argv);
  return calc.Run();
}
#endif
