//=====================================================================-*-C++-*-
// File: $Id: WorkspaceCalculatorConfig.h 581930 2014-02-06 14:20:26Z adye $
//==============================================================================

#ifndef WorkspaceCalculatorConfig_h
#define WorkspaceCalculatorConfig_h

#include "RVersion.h"

#define USE_GLOBITER
#define USE_WildcardList
//#define USE_ProfiledLikelihoodRatioTestStatExt
#define USE_ProfileLikelihoodTestStatEnhanced

#if ROOT_VERSION_CODE>=ROOT_VERSION(5,32,1)
#define ROOT53201
#endif

#if (ROOT_VERSION_CODE==ROOT_VERSION(5,32,0) && ROOT_SVN_REVISION==41754) || \
     ROOT_VERSION_CODE>=ROOT_VERSION(5,33,2)
#define USE_ToyMCImportanceSampler
#define USE_DetailedOutput
#define USE_SetConditionalMLEs
#endif

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
#define USE_StoreFitInfo
#define USE_Analytic
//#define USE_ToyMCImportanceSampler_profile
#endif

#endif
