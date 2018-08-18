/*==============================================================================
$Id: WorkspaceCalculator.cxx 760091 2016-07-06 18:02:47Z adye $

Base class to perform cross-section limits and significance calculations using toys.
This class implements the overall control and setup methods. Derived classes
(LimitCalculator and PValueCalculator) run the toys and plot the results.

Main parameters:

  fileName         workspace file or list of result files
  wsName           workspace name
  modelSBName      ModelConfig object name
  dataName         observed dataset

  calculatorType = 0 Frequentist calculator
                 = 1 Hybrid calculator
                 = 2 Asymptotic calculator
                 = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)

  testStatType   = 0 Simple Likelihood Ratio (LEP)
                 = 1 Ratio of Profile Likelihood (Tevatron)
                 = 2 Profile Likelihood Ratio two-sided
                 = 3 Profile Likelihood Ratio one-sided           (i.e.   0  if mu_hat >= mu)
                 = 4 Maximum Likelihood Estimate (mu_hat)
                 = 5 Profile Likelihood Ratio one-sided discovery (i.e.   0  if mu_hat <= mu)
                 = 6 Profile Likelihood Ratio signed discovery    (i.e. -q0  if mu_hat <  mu)
                 = 7 Profile Likelihood Ratio signed              (i.e. -qmu if mu_hat >  mu)
                 = 8 Number of observed events as test statistic
                 =11 Ratio of Profile Likelihood (Tevatron), subtracting global MLE

  ntoys:           number of toys to use

  plotResult:      plot result of test (TS distribution)
  writeResult:     write result (default is true)

Author: Tim Adye, based on $ROOTSYS/tutorials/roostats/StandardHypoTestInvDemo.C
                  from ROOT trunk revision 41812 (5.32.00-rc1+).

==============================================================================*/

#include "Utilities/WorkspaceCalculatorConfig.h"
#include "WorkspaceCalculator.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <limits>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TError.h"
#include "TString.h"
#include "TPRegexp.h"
#include "TList.h"
#include "TObjArray.h"
#include "TKey.h"
#include "RooRandom.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TText.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TTreeFormula.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStopwatch.h"

#include "Math/MinimizerOptions.h"
#include "Math/IOptions.h"
#include "Math/DistFunc.h"
#include "Minuit2/MnMachinePrecision.h"

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooMsgService.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooRealSumPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooBifurGauss.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"

#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/AsymptoticCalculator.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestResult.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"

#ifdef  USE_GLOBITER
#include "Utilities/GlobIter.h"
#endif

#ifdef  USE_WildcardList
#include "Utilities/WildcardList.h"
#else
struct WildcardList {  // dummy version to keep compiler happy
  WildcardList         (const TString&) {}
  const TString* Match (const TString&) const {return 0;}
  Int_t MatchInd       (const TString&) const {return -1;}
};
#endif

#ifdef USE_ProfiledLikelihoodRatioTestStatExt
#include "cmsfit/interface/ProfiledLikelihoodRatioTestStatExt.h"
#endif

#ifdef USE_ProfileLikelihoodTestStatEnhanced
#include "ProfileLikelihoodTestStatEnhanced.h"
#endif

#ifdef USE_ToyMCImportanceSampler
#include "RooStats/ToyMCImportanceSampler.h"
#endif

#ifdef USE_Analytic
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#endif

#include "Atlas/AtlasStyle.C"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::min;
using std::max;
using std::isnan;

TTime progStart;

const ULong64_t WorkspaceCalculator::default_optimize = (
  WorkspaceCalculator::kReuseNLL            |
  WorkspaceCalculator::kMinuit2             |
  WorkspaceCalculator::kStrategy0           |
  WorkspaceCalculator::kVectorStore         |
  WorkspaceCalculator::kSkipDataTS          |
  WorkspaceCalculator::kUseCondSnapshot     |
  WorkspaceCalculator::kLoadInitialSnapshot |
  WorkspaceCalculator::kPLTSEnhanced
);

// optimize flag characters (-Ox) must correspond to optimize_names and OptimizeBits enums.
static const char* optimize_chars = "RMFSVBIACDsimeru2cLaEgGNfxtWXPbopvOdU";
static const char* const optimize_names[]= {
  "ReuseNLL",
  "UseMultiGen",
  "Minuit2",
  "Strategy0",
  "VectorStore",
  "GenerateBinned",
  "InitialFit",
  "AdjustRanges",
  "SetConst",
  "SkipDataTS",
  "UseCondSnapshot",
  "LoadInitialSnapshot",
  "MinosData",
  "InterpolateErrors",
  "SaveInfo",
  "UseSavedToys",
  "Strategy2",
  "ChainNtuple",
  "AddSLRTS",
  "Analytic",
  "PLTSEnhanced",
  "InitGlobObs",
  "ExpandRanges",
  "DisableNumInt",
  "FixStatError",
  "FixInitially",
  "ReuseAltToys",
  "SaveWorkspace",
  "NoFits",
  "AddPseudoData",
  "BinnedLikelihoodOpt",
  "SetGlobalObservable",
  "FixCache",
  "AltMinos",
  "FastInit",
  "AddErrors",
  "UpdGlobalObservables",
};

TOwnedList::TOwnedList() : TList() { SetOwner(); }
TOwnedList::~TOwnedList()          { Clear(); }
void TOwnedList::Clear(Option_t *option)
{
  if (!option || strcmp(option,"nodelete")!=0)
    for (TIter it(this); TObject* obj= it();) {
//    cout << "Delete "<<obj->ClassName()<<"::"<<obj->GetName()<<(obj->IsOnHeap()?"":" (not on heap)")<<endl;
      delete obj;
    }
  TList::Clear("nodelete");
}

TString WorkspaceCalculator::Join (const std::vector<double>& v, const TString& sep, bool trim)
{
  TString s;
  size_t n=v.size();
  if (trim)
    for (; n>0; n--)
      if (!isnan(v[n-1])) break;
  for (size_t i=0; i<n; i++) {
    if (i>0) s += sep;
    if (!isnan(v[i])) s += Form("%g",v[i]);
  }
  return s;
}

TString WorkspaceCalculator::Join (const TVectorD& v, const TString& sep, bool trim)
{
  TString s;
  size_t n=v.GetNrows();
  if (trim)
    for (; n>0; n--)
      if (!isnan(v[n-1])) break;
  for (size_t i=0; i<n; i++) {
    if (i>0) s += sep;
    if (!isnan(v[i])) s += Form("%g",v[i]);
  }
  return s;
}

TString WorkspaceCalculator::Join (const std::vector<double>& v1, const std::vector<double>& v2,
                                   const TString& sep1, const TString& sep2, bool trim)
{
  TString s;
  size_t n1=v1.size(), n2=v2.size(), n=max(n1,n2);
  if (trim)
    for (; n>0; n--)
      if (!(n>n1 || isnan(v1[n-1])) || !(n>n2 || isnan(v2[n-1]))) break;
  for (size_t i=0; i<n; i++) {
    if (i>0) s += sep2;
    if (isnan(v1[i]) && isnan(v2[i])) continue;
    if (i<n1 && !isnan(v1[i])) s += Form("%g",v1[i]);
    if (i<n2 && isnan(v2[i]) && std::signbit(v2[i])) continue;  // single number from Scan() with ok1=true
    s += sep1;
    if (i<n2 && !isnan(v2[i])) s += Form("%g",v2[i]);
  }
  return s;
}

TString WorkspaceCalculator::Join (const std::vector<int>& v, const TString& sep)
{
  TString s;
  for (size_t i=0, n=v.size(); i<n; i++) {
    if (i>0) s += sep;
    s += Form("%d",v[i]);
  }
  return s;
}

TString WorkspaceCalculator::Join (const std::vector<std::string>& v, const TString& sep, bool trim)
{
  TString s;
  size_t n=v.size();
  if (trim)
    for (; n>0; n--)
      if (v[n-1].size()>0) break;
  for (size_t i=0; i<n; i++) {
    if (i>0) s += sep;
    s += v[i].c_str();
  }
  return s;
}

TString WorkspaceCalculator::Join (const RooAbsCollection& c, const TString& sep, int show,
                                   const TString& sep0, const TString& sep1, const TString& flagc, const TString& flagu)
{
  // show: 0: value, 1: name=value, 2: name=value:min:max, 3:name
  TString s;
  int i= 0;
  for (RooLinkedListIter it= c.iterator(); TObject* o= it.Next();) {
    if ((i++)>0) s += sep;
    if (show>=1) {
      s += o->GetName();
      if (show==3) continue;
      s += sep0;
    }
    if (RooAbsReal* r= dynamic_cast<RooAbsReal*>(o)) {
      if (!isnan(r->getVal())) s += Form("%g",r->getVal());
      if (show>=2) {
        if (r->isConstant()) s += flagc;
        if (RooRealVar* v= dynamic_cast<RooRealVar*>(o)) {
          if (!isnan(v->getMin()) || !isnan(v->getMax())) {
            if (!isnan(v->getVal()) || v->isConstant()) s += sep1;
            if (!isnan(v->getMin())) s += (v->hasMin() ? Form("%g",v->getMin()) : flagu.Data());
            s += sep1;
            if (!isnan(v->getMax())) s += (v->hasMax() ? Form("%g",v->getMax()) : flagu.Data());
          }
        }
      }
    } else
      s += "*";
  }
  return s;
}

const char* WorkspaceCalculator::Scan (const char* a, std::vector<double>& v, char sep, bool append)
{
  if (!append) v.clear();
  if (!*a) return a;
  for (;; a++) {
    if (*a==sep) v.push_back (NaN);
    else         v.push_back (Strtod(a));
    if (*a!=sep) return a;
  }
}

const char* WorkspaceCalculator::Scan (const char* a, std::vector<double>& v1, std::vector<double>& v2,
                                       char sep1, char sep2, bool ok1, bool append)
{
  // With the default sep1=':', sep2=',': specify v1:v2,v1:,:v2, (not just v1, unless ok1=true).
  if (!append) {
    v1.clear();
    v2.clear();
  }
  if (!*a) return a;
  for (size_t n=0;; n++, a++) {
    v1.push_back (NaN);
    v2.push_back (NaN);
    const char* ai= a;
    if (*a==sep2)                   continue;
    if (*a!=sep1 && *a)             v1[n]= Strtod(a);
    if (ok1) {
      if (*a!=sep1 && !isnan(v1[n]))
        v2[n]= copysign(v2[n],-1.0);  // flag single v1 with v2=-NaN (test with std::signbit(v2))
      if (*a==sep2)                 continue;
    } else if (!*a)                 return ai;
    if (*a!=sep1)                   return a;
    a++;
    if (*a!=sep1 && *a!=sep2 && *a) v2[n]= Strtod(a);
    if (*a!=sep2)                   return a;
  }
}

const char* WorkspaceCalculator::Scan (const char* a, RooAbsCollection& vars,
                                       char sep0, char sep1, char sep2, char flagc, char flagn, char inf, bool append)
{
  // With the default sep0='=', sep1=':', sep2=',': specify name=val,name=min:max,name=val:min:max,name=val:min:max:err.
  // Specify, eg., "2.0C" to set constant. Use "-" (or -/+1e30) for lower/upper limit to remove limit.
  // Each name can itself be a comma-separated list, where names beginning with "-" exclude from the rest of that list.
  // The settings are added to the RooAbsCollection vars ('append' is true by default), and applied in SetParameters().
  const char allsep[]= {sep0, sep1, sep2, '\0'};
  if (!append) vars.removeAll();
  flagc= toupper(flagc);
  RooArgList set;
  string excl;
  do {
    const char* b= strpbrk (a, allsep);
    if (!b || *b==sep1) return a;
    string name(a,b);
    RooRealVar* var=0;
    if (!name.empty() && name[0]==flagn) {
      if (excl.empty()) name.erase(0,1);
      else              name[0]= ',';
      excl += name;
      if (*b==sep2) {
        a=b;
        continue;
      }
      if (set.getSize()==0) name= "*";
      else                  var= dynamic_cast<RooRealVar*>(set.at (set.getSize()-1));
    }
    if (!var) {
      var= new RooRealVar (name.c_str(), name.c_str(), NaN);
      var->setRange (NaN, NaN);
      var->setConstant(kFALSE);
      vars.addOwned(*var);
      if (*b==sep2) {
        set.add(*var);
        a=b;
        continue;
      }
      if (!excl.empty()) var->setStringAttribute ("excludeNames", excl.c_str());
    }
    a= b+1;
    Double_t num[4]= {NaN,   NaN,   NaN,   NaN};
    bool   isInf[3]= {false, false, false}, isConst= false;
    int n=0;
    while (n < (isInf[0] ? 2 : 4)) {
      if (n==0 && toupper(*a)==flagc) {
        isConst= true;
        a++;
      } else if (n<3 && *a==inf && (a[1]==sep1 || a[1]==sep2 || !a[1])) {
        isInf[n]= true;
        a++;
      } else if (*a!=sep1 && *a!=sep2 && *a) {
        num[n]= Strtod(a);
        if (n==0 && toupper(*a)==flagc) {
          isConst= true;
          a++;
        }
      }
      n++;
      if (*a!=sep1) break;
      a++;
    }
    if (isInf[0] || (n==2 && !isConst)) {
      var->setRange (num[0], num[1]);
      if (isInf[0]) var->removeMin();
      if (isInf[1]) var->removeMax();
    } else {
      if (n>=1) var->setVal (num[0]);
      if (n>=2) var->setRange (num[1], num[2]);
      if (isInf[1]) var->removeMin();
      if (isInf[2]) var->removeMax();
      if (n>=4) var->setError (num[3]);
      if (isConst) var->setConstant();
    }
    for (RooLinkedListIter it= set.iterator(); RooRealVar* v= dynamic_cast<RooRealVar*>(it.Next());) {
      v->setConstant (var->isConstant());
      v->setVal      (var->getVal());
      v->setRange    (var->getMin(), var->getMax());
      if (var->hasError()) v->setError (var->getError());
      if (!excl.empty()) v->setStringAttribute ("excludeNames", excl.c_str());
    }
    set.removeAll();
    excl= "";
  } while (*a++==sep2);
  return a-1;
}

const char* WorkspaceCalculator::Scan (const char* a, std::vector<int>& v, char sep, int nullVal, bool append)
{
  if (!append) v.clear();
  if (!*a) return a;
  for (;; a++) {
    const char* a0=a;
    int val= Strtol(a);
    if (a==a0) val= nullVal;
    v.push_back (val);
    if (*a!=sep) return a;
  }
}

const char* WorkspaceCalculator::Scan (const char* a, std::vector<std::string>& v, char sep, bool append)
{
  if (!append) v.clear();
  if (!*a) return a;
  for (;; a++) {
    const char* b= strchr (a, sep);
    if (!b) b= a+strlen(a);
    v.push_back (string(a,b));
    if (!*b) return b;
    a= b;
  }
}


size_t WorkspaceCalculator::ReadFile (const char* filename, std::vector<std::string>& v, bool strip, char comment, bool append)
{
  if (!append) v.clear();
  std::ifstream f(filename);
  if (!f) {
    cerr << "Could not read " << filename << endl;
    return 0;
  }
  size_t n=1;  // offset by 1 so 0=fail
  while (f) {
    string s;
    std::getline (f, s);
    if (strip) {
      if (comment!='\0') {
        size_t ic= s.find('#');
        if (ic!=string::npos) s.erase(ic);
      }
      size_t i1= s.find_first_not_of(' ');
      if (i1==string::npos) continue;   // blank line
      size_t i2= s.find_last_not_of(' ');
      if (i1>0 || i2+1<s.size()) s= s.substr(i1,i2-i1+1);
    }
    v.push_back(s);
    n++;
  }
  return n;
}

std::string WorkspaceCalculator::ReadFile (const char* filename)
{
  string contents;
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (!in) {
    cerr << "Could not read " << filename << endl;
    return contents;
  }
  in.seekg(0, std::ios::end);
  contents.resize(in.tellg());
  in.seekg(0, std::ios::beg);
  in.read (&contents[0], contents.size());
  in.close();
  return contents;
}

inline
Bool_t MoveBin (TH1* h, Int_t bin1, Int_t bin2)
{
  Double_t v= h->GetBinContent(bin1);
  if (v==0.0) return kFALSE;
  h->AddBinContent (bin2, v);
  h->SetBinContent (bin1, 0.0);
  if (h->GetSumw2N()) {
    Double_t e1= h->GetBinError(bin1), e2= h->GetBinError(bin2);
    h->SetBinError (bin2, std::sqrt (e1*e1 + e2*e2));
  }
  return kTRUE;
}

Int_t WorkspaceCalculator::IncludeOverflows (TH1* h)
{
  if (!h) return 0;
  Int_t nover=0;
  Double_t entries= h->GetEntries();
  Double_t stats[TH1::kNstat];
  h->GetStats(stats);
  if        (h->GetDimension()==1) {
    Int_t nx= h->GetNbinsX();
    if (MoveBin (h,    0,  1)) nover++;
    if (MoveBin (h, nx+1, nx)) nover++;
  } else if (h->GetDimension()==2) {
    Int_t nx= h->GetNbinsX(), ny= h->GetNbinsY();
    for (Int_t i=0; i<=nx+1; i++) {  // corner bins first get added to x-overflow bins
      if (MoveBin (h, h->GetBin(i,   0), h->GetBin(i, 1))) nover++;
      if (MoveBin (h, h->GetBin(i,ny+1), h->GetBin(i,ny))) nover++;
    }
    for (Int_t i=1; i<=ny; i++) {
      if (MoveBin (h, h->GetBin(0,   i), h->GetBin(1, i))) nover++;
      if (MoveBin (h, h->GetBin(nx+1,i), h->GetBin(nx,i))) nover++;
    }
  } else if (h->GetDimension()==3) {
    Int_t nx= h->GetNbinsX(), ny= h->GetNbinsY(), nz= h->GetNbinsZ();
    for (Int_t i=0; i<=nx+1; i++) {  // edge bins first get added to x,y-overflow bins
      for (Int_t j=0; j<=ny+1; j++) {
        if (MoveBin (h, h->GetBin(i,j,   0), h->GetBin(i,j, 1))) nover++;
        if (MoveBin (h, h->GetBin(i,j,nz+1), h->GetBin(i,j,nz))) nover++;
      }
    }
    for (Int_t i=0; i<=nx+1; i++) {
      for (Int_t j=0; j<=nz+1; j++) {
        if (MoveBin (h, h->GetBin(i,   0,j), h->GetBin(i, 1,j))) nover++;
        if (MoveBin (h, h->GetBin(i,ny+1,j), h->GetBin(i,ny,j))) nover++;
      }
    }
    for (Int_t i=1; i<=ny; i++) {
      for (Int_t j=1; j<=nz; j++) {
        if (MoveBin (h, h->GetBin(   0,i,j), h->GetBin( 1,i,j))) nover++;
        if (MoveBin (h, h->GetBin(ny+1,i,j), h->GetBin(ny,i,j))) nover++;
      }
    }
  }
  if (nover) {
    // restore statistics and entries modified by SetBinContent
    h->SetEntries (entries);
    h->PutStats (stats);
  }
  return nover;
}

Int_t WorkspaceCalculator::RemoveConstantParameters (RooAbsCollection* coll, const RooWorkspace* ws, const char* msg)
{
  // Remove elements from coll that are constant (in the workspace or in coll itself)
  RooArgSet rem;
  for (RooFIter it= coll->fwdIterator(); const RooAbsArg* a= it.next();) {
    if (a->isConstant())
      rem.add(*a);
    else if (ws)
      if (RooAbsArg* aw= ws->arg(a->GetName()))
        if (aw->isConstant())
          rem.add(*a);
  }
  Int_t n= rem.getSize();
  if (n<=0) return n;
  coll->remove (rem);
  if (msg) {
    cout << msg << endl;
    rem.Print("v");
  }
  return n;
}

void WorkspaceCalculator::Init() {
  optimize           = default_optimize;
  plotResult         = 0;     // plot results? 0=no, 1=with captions, 2=without captions
  writeResult        = true;  // write result to a file
  nworkers           = 0;     // set non-zero to use PROOF Light when using toys (for freq or hybrid)
  noSystematics      = false; // force all systematics to be off (i.e. set all nuisance parameters as constat
                              // to their nominal values)
  initialFit         = -1;    // do a first  fit to the model (-1 : default, 0 skip fit, 1 always do fit)
  ntoys              = 1000;  // number of null toys
  ntoysAlt           = -1;    // number of alt toys
  nToysRatio         = -1;    // ratio Ntoys S+B/B. Default is for no alt toys.
  nToysResample      = 50;    // number of toys to use in resampling error calculation
  nToysLimit         = false; // limit number of toys in read sample
  dropNegative       = false; // drop negative test statistic values?
  dropNaN            = true;  // drop NaN test statistic values (can occur with ROOT 5.33.02 ToyMCSampler)?
  negVal             = -0.025;// actually only count values less than this as negative
  dropLarge          = -1.0;  // Drop TS values larger than this
  invMass            = 0.0;   // Only required for final result ntuple
  fitTol             = -1.0;  // Migrad tolerance: target EDM=0.001*fitTol. RooFit default is 1 (different from Minuit default tolerance of 0.01), giving convergence when EDM<0.001.
  fitTol0            = -1.0;  // Initial fit with larger tolerance.
  fitPrec            = -1.0;  // Machine precision for Minuit, relative to Minuit2 default (~8.9e-16)
  OneSidedPositiveMin= -0.1;  // Use this poimin with ProfileLikelihoodTestStat::SetOneSidedPositive (testStatType=5)
  plotFile           = Form("%s.pdf",GetName());
  verbose            = 0;
  doStatsTree        = 1;
  detailedOutput     = 0;     // save detailed output from test statistic (2=withErrorsAndPulls)
  skipPlot           = 0;
  plotMuHatPositive  = false; // Plot NP only for muhat>=0
  bPOIdefault        = NaN;
  sbPOIdefault       = NaN;
  force              = false; // don't exit on consistency error
  initialSnapshot    = "ucmles";  // loaded if optimize|kLoadInitialSnapshot (can be comma-sepaarted list)
  initialSnapshotTS  = "";    // initial snapshot before performing fits in ProfileLikelihoodTestStatEnhanced
  newParameterRanges = 10.0;  // #sigma to use for parameter range adjustment
  parmProf           = -1;    // plot parameter profiles: 0=just data, >0: also for a number of toys
  maxFunctionCalls   = 0;     // RooMinimizer uses 500*N,
                              // raw ROOT::Fit::Fitter/Minuit1/2 uses 200+100*N+5*N*N if 0 (N=num fit variables)
                              // NB. RooMinimizer (unless patched) default is 500*N regardless of the MinimizerOptions default,
                              // so we reset it in ProfileLikelihoodTestStatEnhanced.
  minosSetNum        = 0;     // set number of parameters to determine using MINOS
  minosSetSize       = 0;     // number of paramters to run through MINOS at a time (0=all)
  minosSet           = 0;     // set of parameters to determine using MINOS
  minos_results      = 0;
  skipToys           = 0;     // skip fit for this number of toys
  seed               = 0;     // use seed based on job num, else a unique (UUID) seed each job (see TRandom3::SetSeed)
  fixStatErrorMin    = 0.08;  // fix histogram bin statistics nuisance parameters when error is less than this (and -Of specified)

  statsTree=0;
  poihat=NaN; poierr=0.0;
  doneRead= false;
  testStatType=-1;            // test statistic type
  overrideTestStatType=-1;    // >=0: multiple test stat types in different data files or specified by user
  samplerType=-1;             // -1:no importance sampling, 0:adaptive importance sampling, >0: importance sampling with this number of alt sets
  showTitles         = false;
  showCaptions       = true;
  atlasStyle         = 0;
  haveWeights= 0;
  nStdDevOverlap= 0.5;        // overlap of importance samples
  wsInfo= 0;
  asimovSig          = 0.0;
  varSettings= new RooArgList("settings");
  kept.Add(varSettings);

  wsName             = "combWS";
  modelSBName        = "ModelConfig";
  dataName           = "combData";
  calculatorType     = 0;
  useNumberCounting  = -1;
  jobNum             = -1;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,32,1)
  interpCodeNorm     = 4;
#else
  interpCodeNorm     = 1;
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
  interpCodeShape    = 4;
#else
  interpCodeShape    = 0;
#endif
  altResultName      = false;
  showAtlasLabel     = true;
  runLimit = runToys = runScan = optionalResult = constrainedOnly = false;
  resetGlobsSnapshot = 0;

  initErr            = 0;
  ws                 = 0;

  InitCalculator();
}

void WorkspaceCalculator::InitCalculator() {
  data               = 0;
  sbModel = bModel = nullModel = altModel = 0;
  sbPdf   = bPdf   = nullPdf   = altPdf   = 0;
  nullSnapshot = altSnapshot = 0;
  testStat = testStat2 = 0;
  allParams = nullParams = minosParams = minosSet = 0;
  minos_results      = 0;
  toymcs             = 0;
  calc               = 0;
}

void WorkspaceCalculator::SetDefaults() {
  if (testStatType<0) testStatType= 2;            // test statistic type (2=two-sided)
  if (BkgIsAlt()) {
    nullModelName= "S+B";
     altModelName= "B";
  } else {
    nullModelName= "B";
     altModelName= "S+B";
  }
}

void WorkspaceCalculator::DeleteCalculator()
{
  delete calc;
  delete toymcs;
  delete testStat;
  delete testStat2;
  delete nullSnapshot;
  delete  altSnapshot;
  delete  allParams;
  delete nullParams;
  delete minosParams;
  delete minosSet;
  calcObjs.Clear();
  InitCalculator();
}


int WorkspaceCalculator::SetOptimize (const char* a, const char* opt) {
  if (!a || !*a) return 1;
  char* ae;
  ULong64_t newopt= strtol (a, &ae, 10);
  if (*ae=='\0' && !(a[0]=='2' && a[1]=='\0'))
    optimize= newopt;
  else {
    bool neg= false;
    for (;*a; a++) {
      if      (*a=='-') neg= true;
      else if (*a=='+') neg= false;
      else if (const char* ac= strchr (optimize_chars, *a)) {
        if (neg) optimize &= ~(1LL<<(ac-optimize_chars));
        else     optimize |=  (1LL<<(ac-optimize_chars));
      } else {
        cerr << GetName() << ": invalid option: " << opt << *a << endl; return 0;
      }
    }
  }
  if (optimize & kStrategy2) optimize &= ~kStrategy0;  // for clarity in PrintOptimize message
  return 1;
}

void WorkspaceCalculator::PrintOptimize (ULong64_t opt, const char* first) {
  if (!first && !opt)
    cout << '0';
  else {
    for (int i=0, j=0, l=strlen(optimize_chars); i<l; i++) {
      if (opt & (1LL<<i)) {
        if (first) cout << (j++ ? ", " : first) << optimize_names[i];
        else       cout << optimize_chars[i];
      }
    }
  }
}


int WorkspaceCalculator::Run()
{
   if (initErr) return initErr;
   progStart = gSystem->Now();

   if (verbose<0) gErrorIgnoreLevel= kWarning;
   if (plotResult>=2) {
     showTitles= showCaptions= false;
     atlasStyle= plotResult-1;
   }

   if (seedName.Length()>0) {
     TString name= seedName;
     int i= name.Last(':');
     if (i == kNPOS) i= name.Length();
     TFile* f= TFile::Open (TString(name(0,i)));
     if (!f) return 5;
     if (!ReadObject (Form("rnd/toy_%s", (i+1<name.Length() ? name.Data()+i+1 : "0")),
                      RooRandom::randomGenerator(), f)) return 5;
     delete f;
   } else {
     if (seed < 0) seed= 4357 + jobNum + 53475*int(10.0*invMass+3.5) + 12073*(seed+1);
     RooRandom::randomGenerator()->SetSeed(seed);  // the default, seed=0, uses a unique (UUID) seed each job (see TRandom3::SetSeed)
     if (seed ==0) optimizeTS= "S"+optimizeTS;  // save seeds
   }
   if (verbose>=0)
     cout << "Random number seed = " << RooRandom::randomGenerator()->GetSeed() << endl;  // only the first of 8 for SetSeed(0)

   if (ntoysAlt<0) {
     ntoysAlt= nToysRatio<=0.0 ? 0 :
               nToysRatio==1.0 ? ntoys
                               : TMath::CeilNint(ntoys/nToysRatio);
   }

   // Need initial fit for LEP test stat or adaptive importance sampling, but maybe snapshot is OK
   if (!(optimize & (kLoadInitialSnapshot|kInitialFit)) && (testStatType==0 || samplerType==0))
     optimize |= kLoadInitialSnapshot;

   if (fileName.empty()) {
     cerr << "No file name specified" << endl;
     return 1;
   }

   if (!(optimize & kPLTSEnhanced)) {
     if ((optimize & (kSaveInfo|kUseSavedToys|kMinosData)) || parmProf>=0 || optimizeTS.Length()>0 || initialSnapshotTS.Length()>0)
       cerr << "WARNING: disabling options specific to ProfileLikelihoodTestStatEnhanced" << endl; // also SkipDataTS, which is default
     optimize &= ~(kSaveInfo|kUseSavedToys|kMinosData|kSkipDataTS);
   }

   if (verbose>=0) {
     ShowParms();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,8)
     cout << Form ("ROOT %s (%s@%s, %s on %s)",
                   gROOT->GetVersion(), gROOT->GetGitBranch(), gROOT->GetGitCommit(),
                   gROOT->GetGitDate(), gSystem->GetBuildArch()) << endl;
#else
     cout << Form ("ROOT %s (%s@%d, %s on %s)",
                   gROOT->GetVersion(), gROOT->GetSvnBranch(), gROOT->GetSvnRevision(),
                   gROOT->GetSvnDate(), gSystem->GetBuildArch()) << endl;
#endif
   }

   if (optimize & kVectorStore) RooAbsData::setDefaultStorageType (RooAbsData::Vector);
   else                         RooAbsData::setDefaultStorageType (RooAbsData::Tree);

   bool singleFile=false, haveResult= false;
   if (fileName.size()==1 && TString(fileName[0]).First("[]*?") == kNPOS) {  // just one file, no wildcards
     singleFile= true;
     TFile* file = TFile::Open (fileName[0].c_str());
     if (!file) {
       cerr << "Could not open file " << fileName[0] << endl;
       return 2;
     }
     file->GetObject (wsName, ws);
     gROOT->cd();
     if (ws) {
       kept.Add (file);
       wsfile= fileName[0];
       SetDefaults();
       Info("Run","Use workspace '%s' in %s",wsName.Data(),fileName[0].c_str());
       int ok= RunOverWorkspace();
       if (!saved.FindObject(ws)) kept.Add(ws);
       if (!ok) {
         cerr << "Error running the calculator - Exit" << endl;
         return 3;
       }
       haveResult= true;
     } else
       delete file;
   }


   if (!haveResult) {
     if (!ReadResult (fileName)) {
       if (singleFile) Error("Run","No workspace '%s' in file %s",wsName.Data(),fileName[0].c_str());
       return 4;
     }
     doneRead= true;
     SetDefaults();
     if (poimax.empty()) poimax.assign(1,10.0);
     if (poimin.empty()) {
       if      (testStatType==5) poimin.assign(1,OneSidedPositiveMin);
       else if (testStatType==6) poimin.assign(1,-poimax[0]);
       else                      poimin.assign(1,0.0);
     }
     if (nullPOI.empty() && altPOI.empty()) {
       if (BkgIsAlt()) {
         nullPOI      = sbPOI;
          altPOI      =  bPOI;
       } else {
         nullPOI      =  bPOI;
          altPOI      = sbPOI;
       }
     }
   }
   std::cout << "XXX 2nd end of RUN" << std::endl;

   AdjustResults(!haveResult);
   AnalyseResults();

   // Deleting a large workspace is much slower than destroying it on exit
   // This does mean a memory leak if we delete WorkspaceCalculator without exiting.
   kept.Remove(ws);
   saved.Remove(ws);
   std::cout << "XXX end of RUN" << std::endl;

   if (verbose>=0) PrintResourcesUsed(progStart);
   return 0;
}


void WorkspaceCalculator::AdjustResults (bool resultWasRead)
{
  if (resultWasRead && nToysLimit) LimitSamples();
  DropBadSamples();
  ApplyCuts();
}


void WorkspaceCalculator::AnalyseResults()
{
  if (verbose>=0) PrintResult();
  GetLimits();
  if (plotResult) PlotResult();
  if (writeResult) WriteResultFile();
}

void
WorkspaceCalculator::PlotResult()
{
  if (atlasStyle) {
    SetAtlasStyle(false);
    gStyle->SetHistLineColor(kBlue+2);
  }
  gStyle->SetMarkerColor(gStyle->GetHistLineColor());
  if (!atlasStyle) gStyle->SetTitleYOffset(1.2);
  gStyle->SetOptStat (0);
  gStyle->SetOptTitle (showTitles ? 1 : 0);
  if (!showTitles)                  gStyle->SetPadTopMargin   (0.03);
  if (!showCaptions && !atlasStyle) gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.14);
  cout << "Plot results with " << gStyle->GetTitle() << endl;

  TCanvas* canvas = new TCanvas(GetName(),GetName(),800,(atlasStyle==1?600:800));
  canvas->SetFillStyle(0);
  canvas->Print (Form("%s[",plotFile.Data()));

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);

  PlotResult(canvas);

  TH1::AddDirectory (oldstat);

  canvas->Print (Form("%s]",plotFile.Data()));
  delete canvas;
}

int WorkspaceCalculator::WriteResultFile()
{
  if (!(writeResult && !resultFileName.empty() && resultFileName!="-" && saved.First())) return 0;

  if (!force) {
    for (size_t j=0; j<fileName.size(); j++) {
      if (resultFileName!=fileName[j]) continue;
      cerr << "Do not overwrite input file "<<resultFileName<<" - skip writing result file"<<endl;
      return 0;
    }
  }

  if (resetGlobsSnapshot && ws) {
    cout << "Reset initial global observables from before Asimov generation" << endl;
    ws->loadSnapshot(resetGlobsSnapshot);
  }

  Info("Run","Write result to file %s",resultFileName.c_str());
  TFile* fileOut = TFile::Open (resultFileName.c_str(), "recreate");
  if (!fileOut) return 0;
  for (TIter it= &saved; TObject* obj= it();) {
    TString objName, objPath;
    // If TObject has default GetName (returning class name), use option if set.
    // Used for TMatrixDSym saved in ProfileLikelihoodTestStatEnhanced
    if (obj->GetName() == obj->IsA()->GetName() && *it.GetOption()) objPath= it.GetOption();
    else                                                            objPath= obj->GetName();
    TDirectory* d;
    TString dirName= gSystem->DirName(objPath);
    if (dirName == ".") {
      d= fileOut;
      objName= objPath;
    } else {
      d= fileOut->GetDirectory(dirName);
      if (!d) d= fileOut->mkdir(dirName);
      objName= gSystem->BaseName(objPath);
    }
    TTree* tree= dynamic_cast<TTree*>(obj);
    if (tree && tree==statsTree && cut.Length()>0) {
      TString ntcut;
      ntcut.Form("(%s) || is_data",cut.Data());
      Info("Run","Write tree '%s' to file %s, applying cut '%s'",
           tree->GetName(), resultFileName.c_str(), cut.Data());
      d->cd();
      if (TTree* newTree= tree->CopyTree(ntcut)) {
        Long64_t nent= tree->GetEntries(), nsel= newTree->GetEntries();
        if (verbose>=0) cout << nent-nsel << " / " << nent << " entries fail cut" << endl;
        d->WriteTObject (newTree, objName);
      }
      fileOut->cd();
    } else if (TChain* chain= dynamic_cast<TChain*>(obj)) {
      Info("Run","Merge chain '%s' of %d files into file %s",
           chain->GetName(), chain->GetListOfFiles()->GetEntries(), resultFileName.c_str());
      if (verbose>=3) chain->ls();
      d->cd();
      chain->Merge(fileOut,0,"keep");  // with "fast,SortBasketsByEntry" it's faster, but uses too much memory
      fileOut->cd();
    } else {
      if (verbose>=1) {
        cout << "Write "<<obj->ClassName()<<"::"<<objPath<<endl;
        if (verbose>=3) obj->Print("v");
      }
      d->WriteTObject (obj, objName);
    }
  }
  delete fileOut;
  Info("Run","Write OK");
  return 1;
}


const RooArgSet* WorkspaceCalculator::SetValueInSnapshot (RooStats::ModelConfig* mc, const std::vector<double>& val)
{
  RooArgSet* poiSet= (RooArgSet*)mc->GetParametersOfInterest();
  if (!poiSet) return 0;
  RooArgSet* old= 0;
  if (!val.empty()) {
    old= dynamic_cast<RooArgSet*>(poiSet->snapshot());
    size_t i= 0;
    for (RooLinkedListIter it = poiSet->iterator(); TObject* o = it.Next(); i++) {
      if (i >= val.size()) break;
      if (RooRealVar* poi= dynamic_cast<RooRealVar*>(o)) {
        if (val[i] < poi->getMin()) poi->setMin(val[i]);
        if (val[i] > poi->getMax()) poi->setMax(val[i]);
        poi->setVal(val[i]);
      }
    }
  }
  mc->SetSnapshot (*poiSet);
  const RooArgSet* snap= mc->GetSnapshot();
  if (old) {
    *poiSet= *old;
    delete old;
  }
  return snap;
}

Bool_t WorkspaceCalculator::LoadSnapshot (RooWorkspace* w, const TString& snapName, const char* msgfmt)
{
  if (snapName.Length()<=0) return kFALSE;
  Bool_t loaded= kFALSE;
  for (int i=0,j=0; j!=kNPOS; i=j+1) {
    j= snapName.Index(",",i);
    TString name= snapName (i, j!=kNPOS ? j-i : snapName.Length());
    if (!w->loadSnapshot(name)) continue;
    loaded= kTRUE;
    if (msgfmt) ::Info("WorkspaceCalculator::LoadSnapshot",msgfmt,name.Data());
  }
  return loaded;
}


RooArgSet* WorkspaceCalculator::GetSnapshot (RooWorkspace* w, const TString& snapName, const char* msgfmt,
                                             const RooArgSet* vars, const char* resetSnap /* = "initialNuisAndPoi" */)
{
  // Get snapshot(s) from workspace, saving values in a copy of vars.
  // For ROOT 5.34.05 and before, the original workspace state can only be recovered by loading initialNuisAndPoi.
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,6)
  if (resetSnap) {}  // prevent unused parameter warning
  if (snapName.Length()<=0) return 0;
  RooArgSet* snap=0;
  for (int i=0,j=0; j!=kNPOS; i=j+1) {
    j= snapName.Index(",",i);
    TString name= snapName (i, j!=kNPOS ? j-i : snapName.Length());
    const RooArgSet* s= w->getSnapshot(name);
    if (!s) continue;
    if (msgfmt) ::Info("WorkspaceCalculator::GetSnapshot",msgfmt,name.Data());
    if (!snap) snap= dynamic_cast<RooArgSet*>(vars->snapshot());
    *snap= *s;
  }
#else
  if (!LoadSnapshot (w, snapName, msgfmt)) return 0;
  RooArgSet* snap= dynamic_cast<RooArgSet*>(vars->snapshot());
  if (resetSnap) assert (w->loadSnapshot(resetSnap));
#endif
  return snap;
}


void WorkspaceCalculator::ModifyInterpolationForAll(RooWorkspace* w, int codeNorm, int codeShape){
  // Modify the interpolation used in HistFactory models.
  // Codes for interpolation
  //   code = 0: piece-wise linear
  //   code = 1: pice-wise log
  //   code = 2: parabolic interp with linear extrap
  //   code = 3: parabolic version of log-normal
  //   code = 4: polynomial interpolation, log extrapolation
  // currently same codes for shape interpolation, but maybe not forever, so we separate them.
  //
  // From $ROOTSYS/tutorials/histfactory/ModifyInterpolation.C,
  // modified to enable both normalization and shape interpolations together.
  int nn=0, ns=0;
  RooArgSet funcs = w->allFunctions();
  for (RooLinkedListIter it = funcs.iterator(); TObject* tempObj=it.Next();) {
    if (codeNorm>=0) {
      RooStats::HistFactory::FlexibleInterpVar* flex = dynamic_cast<RooStats::HistFactory::FlexibleInterpVar*>(tempObj);
      if(flex){
        flex->setAllInterpCodes(codeNorm);
        if (verbose>=1) flex->printAllInterpCodes();
        nn++;
      }
    }
    if (codeShape>=0) {
      PiecewiseInterpolation* piece = dynamic_cast<PiecewiseInterpolation*>(tempObj);
      if(piece){
        piece->setAllInterpCodes(codeShape);
        if (verbose>=1) piece->printAllInterpCodes();
        ns++;
      }
    }
  }
  if (verbose==0) {
    if (nn) Info("ModifyInterpolationForAll","Set normalization interpolation code %d for %d functions", codeNorm,  nn);
    if (ns) Info("ModifyInterpolationForAll","Set shape         interpolation code %d for %d functions", codeShape, ns);
  }
}


RooDataSet* WorkspaceCalculator::SetDatasetBinning (const RooAbsPdf* pdfIn, const RooAbsData* data,
                                                    Int_t setNbins, const char* generateBinnedTag,
                                                    const char*   binnedCategories,
                                                    const char* unbinnedCategories, int verbose,
                                                    const char* weightVarName)
{
  WildcardList   binnedCategoriesWild (  binnedCategories ?   binnedCategories : "");
  WildcardList unbinnedCategoriesWild (unbinnedCategories ? unbinnedCategories : "");
  TOwnedList localObjs;
  const RooSimultaneous* pdf= dynamic_cast<const RooSimultaneous*>(pdfIn);
  if (!pdf) return 0;
  const RooAbsCategoryLValue& cat= pdf->indexCat();
  TList* dataList= data->split(cat,true);
  if (!dataList) return 0;
  localObjs.Add(dataList);
  std::map<string,RooDataSet*> dataMap;
  RooRealVar weightVar (weightVarName, "", 1.0);
  int dsBinned= 0;
  for (TIter nextds= dataList; RooAbsData* datai= dynamic_cast<RooAbsData*>(nextds());) {
    localObjs.Add(datai);  // make sure we delete them
    // Make a copy just so we can match the weight variable names
    TString copyName= Form("%s_unbinned",datai->GetName());
    RooDataSet* copyData= new RooDataSet (copyName, copyName, RooArgSet(*datai->get(),&weightVar), weightVar.GetName());
    copyData->append (*(RooDataSet*)datai);   // also works for RooDataHist
    dataMap[datai->GetName()]= copyData;
    localObjs.Add(copyData);
    int isBinned= -1;
    const char* dataType= "";
    if (RooAbsPdf* pdfi= pdf->getPdf(datai->GetName()))
      if (RooArgSet* obs= pdfi->getObservables(*datai)) {
        isBinned= pdfi->isBinnedDistribution(*obs) ? 1 : 0;
        dataType= (isBinned ? " binned" : " unbinned");
        delete obs;
      }
    if (verbose>=1) cout << Form("Category %s%s dataset has %d/%g entries",
                                 datai->GetName(), dataType, datai->numEntries(), datai->sumEntries()) << endl;
    if (verbose>=2) {
      datai->get(0);
      datai->Print("v");
    }
    if (datai->numEntries() <= 0) continue;
    bool rebin= false, genBinOnly= false;
    if        (isBinned==1 && unbinnedCategoriesWild.Match(datai->GetName())) {
      if (datai->numEntries() <= datai->sumEntries())    continue;
      rebin= true;
    } else if (isBinned==0 &&   binnedCategoriesWild.Match(datai->GetName())) {
      if (setNbins <= 0)                                 continue;
      genBinOnly= (datai->numEntries() <= setNbins);
      if (genBinOnly && datai->sumEntries() <= setNbins) continue;
    } else
      continue;
    RooAbsPdf* pdfi= pdf->getPdf(datai->GetName());
    if (!pdfi) continue;
    RooArgSet* obs= pdfi->getObservables(*datai);
    if (!rebin && generateBinnedTag && *generateBinnedTag)
      pdfi->setAttribute(generateBinnedTag);
    if (verbose>=0)
      cout << Form("%s binning on%s dataset %s with %d/%g entries, PDF %s, variables",
                   (rebin?"Prune":genBinOnly?"Generate":"Set"),
                   dataType, datai->GetName(), datai->numEntries(), datai->sumEntries(), pdfi->GetName());
    for (RooLinkedListIter it = obs->iterator(); TObject* o= it.Next();) {
      RooRealVar* v= dynamic_cast<RooRealVar*>(o);
      if (!v) continue;
      if (verbose>=0) cout << " " << v->GetName();
      if (!rebin) v->setBinning (RooBinning (setNbins, v->getMin(), v->getMax()));
    }
    if (verbose>=0) cout << endl;
    if (genBinOnly) {
      delete obs;
      continue;
    }
    if (obs->getSize()!=1 && !rebin) {
      cerr << "PDF "<<pdfi->GetName()<<" has "<<obs->getSize()<<" observables - can't create binned dataset" <<endl;
      delete obs;
      continue;
    }
    RooRealVar* obsVar= dynamic_cast<RooRealVar*>(obs->first());

    TString  newName= Form ("%s_binned", datai->GetName());
    RooDataSet* dataiNew= new RooDataSet (newName, newName, RooArgSet(*obs,weightVar), weightVar.GetName());
    localObjs.Add(dataiNew);
    if (rebin) {
      if (verbose>=3) {
        const RooAbsBinning& binning= obsVar->getBinning();
        cout << Form("Observable %s bins (%d,%g,%g) contain %d/%g entries:",
                     obsVar->GetName(), binning.numBins(), binning.lowBound(), binning.highBound(),
                     datai->numEntries(), datai->sumEntries());
      }
      for (int j=0, n=datai->numEntries(); j<n; j++) {
        const RooArgSet* row = datai->get(j);
        Double_t w= datai->weight();
        if (w == 0) continue;
        if (verbose>=3) cout <<" "<<dynamic_cast<RooRealVar*>(row->find(*obsVar))->getVal()<<":"<<w;
        dataiNew->add (*row, w);
      }
      if (verbose>=3) cout << endl;
      if (dataiNew->numEntries()==0)
        dataiNew->add (*datai->get(0), 0.0);  // keep at least one entry (is this necessary?)
      if (dataiNew->numEntries() >= datai->numEntries()) {
        delete obs;
        continue;
      }
      if (verbose>=1)
        cout << Form("Rebin dataset from %d/%g to %d/%g in %s", datai->numEntries(), datai->sumEntries(), dataiNew->numEntries(), dataiNew->sumEntries(), dataiNew->GetName()) << endl;
    } else {
      TString histName= Form ("%s_hist",   datai->GetName());
      TH1* hist= datai->createHistogram (histName, *obsVar);
      if (verbose>=1)
        cout << Form("Created histogram %s(%d,%g,%g)",
                     hist->GetName(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
      for (int j=1, n=hist->GetNbinsX(); j<=n; ++j) {
        Double_t x= hist->GetXaxis()->GetBinCenter(j), y= hist->GetBinContent(j);
        if (y==0.0) continue;
        if (verbose>=3) cout <<" "<<x<<":"<<y;
        obsVar->setVal (x);
        dataiNew->add (*obs, y);
      }
      if (dataiNew->numEntries()==0) {  // keep at least one entry (is this necessary?)
        obsVar->setVal (hist->GetXaxis()->GetBinCenter(1));
        dataiNew->add (*obs, 0.0);
      }
      if (verbose>=1)
        cout << Form(" in dataset %s with %d/%g entries",
                     dataiNew->GetName(), dataiNew->numEntries(), dataiNew->sumEntries()) << endl;
      delete hist;
    }
    delete obs;
    dataMap[datai->GetName()]= dataiNew;  // replace datai
    dsBinned++;
  }

  if (dsBinned<=0) return 0;

  RooCategory* catVar= dynamic_cast<RooCategory*>(data->get()->find(cat));
  if (!catVar) {
    cerr << "Dataset "<<data->GetName()<<" does not have an index variable '"
         << cat.GetName() << "' - can't create binned dataset" << endl;
    return 0;
  }

  RooArgSet newObs, *allObs= pdf->getObservables(*data);
  for (RooFIter it= data->get()->fwdIterator(); const RooAbsArg* dsvar= it.next();)
    if (RooAbsArg* v= allObs->find(*dsvar)) newObs.add(*v);
  delete allObs;
  TString name=  Form("%s_binnedgg",data->GetName());
  RooDataSet* newData= new RooDataSet (name, name,
                                       RooArgSet (newObs, weightVar),
                                       RooFit::Index(*catVar),
                                       RooFit::Import(dataMap),
                                       RooFit::WeightVar(weightVar));

  if (verbose==0)
    cout << Form("Replace dataset %s (%d/%g entries) with dataset %s (%d/%g entries)",
                 data->GetName(),    data->numEntries(),    data->sumEntries(),
                 newData->GetName(), newData->numEntries(), newData->sumEntries()) << endl;
  else if (verbose>=1) {
    data->get(0);
    newData->get(0);
    cout << Form("=== Replace dataset %s (%d/%g entries):",
                 data->GetName(),    data->numEntries(),    data->sumEntries()) << endl;
    data->Print("v");
    cout << Form("=== with    dataset %s (%d/%g entries):",
                 newData->GetName(), newData->numEntries(), newData->sumEntries()) << endl;
    newData->Print("v");

    if (TList* newDataList= newData->split(cat,true)) {
      for (TIter it= newDataList; RooDataSet* datai= dynamic_cast<RooDataSet*>(it());) {
        cout << Form("New dataset for category %s has %d/%g entries",
                     datai->GetName(), datai->numEntries(), datai->sumEntries()) << endl;
        if (datai->numEntries()!=0) {
          RooArgSet* obs= pdf->getPdf(datai->GetName())->getObservables(*datai);
          RooAbsCollection* dsobs= datai->get(0)->selectCommon(*obs);
          obs->Print("v");
          delete dsobs;
          delete obs;
        }
        delete datai;
      }
      delete newDataList;
    }
  }

  return newData;
}


int WorkspaceCalculator::AddPseudoData (RooAbsData* data, const RooAbsPdf* pdfIn,
                                        const char* categories, const char* binning, int verbose, Double_t eps)
{
  vector<string> vbinning;
  vector< vector<Double_t> > bins;
  if (*Scan (binning, vbinning)) {
    cerr << "WorkspaceCalculator::AddPseudoData: bad binning spec: " << binning << endl;
    return 0;
  }
  for (size_t i=0, n=vbinning.size(); i<n; i++) {
    vector<Double_t> vb;
    if (*Scan (vbinning[i].c_str(), vb, ':') || vb.size()!=3 || vb[0]>=vb[1] || vb[2]<=0.0) {
      cerr << "WorkspaceCalculator::AddPseudoData: bad binning spec " << vbinning[i] << " in " << binning << endl;
      return 0;
    }
    bins.push_back(vb);
  }
  if (bins.empty()) return 0;

  Int_t nument= data->numEntries();
  Double_t sument= data->sumEntries();
  WildcardList categoriesWild (categories);
  const RooSimultaneous* pdf= dynamic_cast<const RooSimultaneous*>(pdfIn);
  if (!pdf) return 0;
  int totcount=0, kilroy=0;
  RooArgSet* obs= pdf->getObservables(*data);
  RooAbsCategoryLValue& cat= const_cast<RooAbsCategoryLValue&>(pdf->indexCat());
  TIterator* iter= cat.typeIterator();
  while (const RooCatType* ct= dynamic_cast<const RooCatType*>(iter->Next())) {
    const char* c= ct->GetName();
    if (verbose>=1) cout << "Category " << c << endl;
    if (!categoriesWild.Match(c)) continue;
    RooAbsPdf* pdfi= pdf->getPdf(c);
    RooArgSet* obsi= pdfi->getObservables(data);
    const RooAbsArg* obsv= obsi->first();
    delete obsi;
    if (!obsv) continue;
    RooRealVar* v= dynamic_cast<RooRealVar*>(obs->find(obsv->GetName()));
    if (!v) continue;
    cat.setLabel(c);
    int count= 0;
    if (!kilroy && verbose>=1) cout << "Add pseudo-data points:";
    for (size_t i=0, n=bins.size(); i<n; i++) {
      vector<Double_t>& b= bins[i];
      Double_t xlo=b[0], xhi=b[1], xstep=b[2]; // xlo<=x<xhi
      for (int j=0, nb= int((xhi-xlo)/xstep+0.9999); j<nb; j++, count++) {
        Double_t xv= xlo+xstep*j;
        if (!kilroy && verbose>=1) cout << ' ' << xv;
        v->setVal(xv);
        data->add(*obs,eps);
      }
    }
    if (!kilroy && verbose>=1) cout << endl;
    kilroy++;
    if (verbose>=0)
      cout << "Added " << count << " pseudo points to observable " << v->GetName()
           << " in category " << c << endl;
    totcount += count;
  }
  delete iter;
  delete obs;
  if (totcount) {
    if (verbose>=0)
      cout << Form("Dataset %s updated from %d/%g entries to %d/%g entries",
                   data->GetName(), nument, sument, data->numEntries(), data->sumEntries()) << endl;
    if (verbose>=1) {
      cout << "Updated dataset:";
      data->Print("v");
    }
    return totcount;
  }
  return 0;
}



RooAbsArg* WorkspaceCalculator::client (const RooAbsArg* a, const TClass* cls1, const TClass* cls2)
{
  // Return unique client that matches either class. Specify two classes to return 0 if both are clients.
  RooArgList clients;
  for (RooFIter ic= a->valueClientMIterator(); RooAbsArg* c= ic.next();) {
    if (!(cls1||cls2) || c->IsA()==cls1 || c->IsA()==cls2) clients.add(*c);
  }
  if (clients.getSize()!=1) return 0;
  return clients.first();
}

RooArgSet WorkspaceCalculator::servers (const RooAbsArg* a)
{
  RooArgSet pargs;
  for (RooFIter is= a->serverMIterator(); RooAbsArg* s= is.next();)
    pargs.add(*s);
  return pargs;
}

bool WorkspaceCalculator::CheckGaussian (RooRealVar* v, const RooArgSet* gobs, bool verbose)
{
  RooGaussian* gaussian= dynamic_cast<RooGaussian*>(client (v, RooGaussian::Class(), RooPoisson::Class()));
  if (!gaussian) return false;
  RooArgSet pargs= servers(gaussian);
  if (pargs.getSize()!=3) return false;
  pargs.remove(*v);
  if (pargs.getSize()!=2) return false;
  pargs.remove(*gobs,kFALSE);
  if (pargs.getSize()!=1) return false;
  const RooAbsReal* sig= dynamic_cast<const RooAbsReal*>(pargs.first());
  if (!sig) return false;
  Double_t sigErr= sig->getVal();
  if (v->getVal()-sigErr < v->getMin()) return false;
  if (v->getVal()+sigErr > v->getMax()) return false;
  v->setError(sigErr);
  if (verbose) {
    cout << "Set error ";
    v->Print("i");
    cout << ": ";
    gaussian->Print();
  }
  return true;
}

bool WorkspaceCalculator::CheckPoisson (RooRealVar* v, bool verbose)
{
  RooPoisson* poisson= dynamic_cast<RooPoisson*>(client (v, RooPoisson::Class(), RooGaussian::Class()));
  if (!poisson) return false;
  RooArgSet pargs= servers(poisson);
  if (pargs.getSize()!=2) return false;
  pargs.remove(*v);
  if (pargs.getSize()!=1) return false;
  const RooAbsReal* nom= dynamic_cast<const RooAbsReal*>(pargs.first());
  if (!nom) return false;
  Double_t nomErr= sqrt(nom->getVal());
  if (v->getVal()-nomErr < v->getMin()) return false;
  if (v->getVal()+nomErr > v->getMax()) return false;
  v->setError(nomErr);
  if (verbose) {
    cout << "Set error ";
    v->Print("i");
    cout << ": nominal ";
    nom->Print();
  }
  return true;
}

bool WorkspaceCalculator::CheckGamma (RooRealVar* v, bool verbose)
{
  RooProduct* prod= dynamic_cast<RooProduct*>(client (v, RooProduct::Class()));
  if (!prod) return false;
  RooArgSet comp= prod->components();
  if (comp.getSize()!=2) return false;
  comp.remove(*v);
  if (comp.getSize()!=1) return false;
  const RooAbsReal* nom= dynamic_cast<const RooAbsReal*>(comp.first());
  if (!nom) return false;
  Double_t nomErr= sqrt(1.0/nom->getVal());
  if (v->getVal()-nomErr < v->getMin()) return false;
  if (v->getVal()+nomErr > v->getMax()) return false;
  v->setError(nomErr);
  if (verbose) {
    cout << "Set error ";
    v->Print("i");
    cout << ": nominal ";
    nom->Print();
  }
  return true;
}


void WorkspaceCalculator::AddErrors (RooStats::ModelConfig* mc, bool verbose)
{
  const RooArgSet* nps=  mc->GetNuisanceParameters();
  const RooArgSet* gobs= mc->GetGlobalObservables();
  TPRegexp CMS_stat("^CMS_.*_stat");
  for (RooFIter it= nps->fwdIterator(); RooAbsArg* a= it.next();) {
    RooRealVar* v= dynamic_cast<RooRealVar*>(a);
    if (!v || v->isConstant() || v->hasError(kFALSE)) continue;
    TString name= v->GetName();
    if (name.Contains("gamma_stat_") && v->getVal()==1.0 && CheckGamma (v, verbose)) continue;
    if (v->getVal()==0.0 && abs(v->getMin()+v->getMax())<1e-9) {
      if ((name.BeginsWith("CMS_") || name.BeginsWith("interf_ggH")) && CheckGaussian (v, gobs, verbose)) continue;
      if (v->getMax()>2.0) {
        v->setError(1.0);
        if (verbose) {
          cout << "Set error ";
          v->Print();
        }
        continue;
      }
    }
    if (name.Contains(CMS_stat) && CheckPoisson (v, verbose)) continue;
  }
}


RooArgSet* WorkspaceCalculator::getAllConstraints (const RooAbsPdf& pdf, const RooArgSet& observables, const RooArgSet* np, int verbose)
{
  // Like RooAbsPdf::getAllConstraints, finds and collects all constraints terms of all component p.d.f.s
  // and returns a RooArgSet with all those terms, except that it expands constraints that depend on other constraints.
  // Also can specify np=0 to return all constraints (does not update np set).
  RooArgSet* constrainedParams;
  if (np)
    constrainedParams= new RooArgSet(*np);
  else {
    constrainedParams= pdf.getParameters(observables,kTRUE);
  }
  RooArgSet* found= new RooArgSet (Form("%s_constraints",pdf.GetName()));
  RooArgSet* c1= pdf.getAllConstraints(observables,*constrainedParams,kFALSE);
  if (verbose>=2 && c1->getSize())
    cout << pdf.ClassName() << "::" << pdf.GetName() << " has constraint terms " << c1->contentsString() << endl;
  for (RooFIter it= c1->fwdIterator(); RooAbsArg* a= it.next();) {
    if (const RooAbsPdf* p= dynamic_cast<const RooAbsPdf*>(a)) {
      RooArgSet* c2= getAllConstraints(*p,observables,constrainedParams,verbose);
      found->add (c2->getSize()>0 ? *c2 : *c1);
      delete c2;
    }
  }
  delete c1;
  delete constrainedParams;
  return found;
}


RooArgSet* WorkspaceCalculator::findGlobalObservable (const RooAbsPdf* pdf, const RooArgSet* observables, const RooArgSet* np, const RooArgSet* globalObservables,
                                                      const RooArgSet* constraints, int verbose)
{
  // Find global observables associated with the specified NPs.
  // If specified, the search is limited to the variables listed in globalObservables.
  // If np is not specified (0), then use all non-constant parameters to guess.
  // If calling multiple times for large workspaces (eg HWW couplings), specify constraints list calculated once with getAllConstraints.
  RooArgSet *vars=0, *fvars=0, *cvars=0, *cons=0;
  if (!np) {
    vars= pdf->getParameters(*observables,kTRUE);
    np= fvars= dynamic_cast<RooArgSet*>(vars->selectByAttrib("Constant",kFALSE));
    if (globalObservables) fvars->remove(*globalObservables,kTRUE);
  }
  if (!constraints) {
    constraints= cons= getAllConstraints(*pdf,*observables,np,verbose);
    if (verbose>=1 && constraints->getSize() > np->getSize())
      cout << "More constraint terms "<<constraints->contentsString()<<" than NPs "<<np->getSize()<<endl;
  }
  if (!globalObservables) {
    if (!vars) vars= pdf->getParameters(*observables,kTRUE);
    globalObservables= cvars= dynamic_cast<RooArgSet*>(vars->selectByAttrib("Constant",kTRUE));
    cvars->remove(*np,kTRUE);
  }
  RooArgSet* found= new RooArgSet(Form("%s_globalObservables",pdf->GetName()));
  for (RooFIter it= constraints->fwdIterator(); RooAbsArg* p= it.next();) {
    if (!cons && !p->dependsOn(*np)) continue;  // if we used a generic list of constraints (passed arg), then check that this term depends on our NPs
    RooArgSet* o= p->getObservables(*globalObservables);
    if (!cvars) o->remove(*np,kTRUE);           // if we used a generic list of global observables (passed arg), then we might have one of our NPs in there
    if (o->getSize()>1) {
      if (verbose>=1) cout << "Constraint "<<p->ClassName()<<"::"<<p->GetName()<<" has more than one potential global observable: "<<o->contentsString()<<" - take first"<<endl;
      found->add(*o->first());
    } else
      found->add(*o);
    delete o;
  }
  delete cons;
  delete cvars;
  delete fvars;
  delete  vars;
  return found;
}


RooArgSet* WorkspaceCalculator::SelectConstrained (const RooArgSet* parms,
                                                   const RooAbsPdf* pdf, const RooArgSet* globalObservables,
                                                   bool constrained, int verbose)
{
  const RooArgSet* allPdfs= pdf->getComponents();
  RooArgSet cons;
  for (RooLinkedListIter it= allPdfs->iterator(); const TObject* o= it.Next();) {
    if (o->IsA()!=RooGaussian::Class() && o->IsA()!=RooPoisson::Class() && o->IsA()!=RooBifurGauss::Class()) continue;
    const RooAbsPdf* pdfi= dynamic_cast<const RooAbsPdf*>(o);
    if (!pdfi) continue;
    bool isconstrained= false;
    RooArgSet args;
    for (RooFIter jt= pdfi->serverMIterator(); const RooAbsArg* a= jt.next();) {
      if (const RooRealVar* v= dynamic_cast<const RooRealVar*>(a)) {
        if      (globalObservables->find(*v)) isconstrained= true;
        else if (parms->find(*v))             args.add(*v);
      }
    }
    if (verbose>=1)
      cout << (isconstrained ? "constraint " : "           ") << pdfi->ClassName() << "::" << pdfi->GetName()
           << " -> " << Join(args,",",3) << endl;
    cons.add (args, kTRUE);
  }
  delete allPdfs;

  // For unconstrained parms, make sure they are still dependants of the pdf
  const RooArgSet* obs= constrained ? 0 : pdf->getObservables (parms);

  // Go through again, so we keep original order of parms
  TString name= parms->GetName();
  name += (constrained ? "_constrained" : "_unconstrained");
  RooArgSet* p= new RooArgSet (name);
  for (RooFIter it= parms->fwdIterator(); RooAbsArg* v= it.next();) {
    bool isconstrained= cons.find(*v);
    if (constrained ? isconstrained : (!isconstrained && obs->find(*v))) p->add(*v);
  }

  delete obs;
  return p;
}

int WorkspaceCalculator::NominalValue (const RooWorkspace* ws, const TString& name, Double_t& nomVal, Double_t& nomErr)
{
  nomVal= nomErr= NaN;
  if (name.Contains("gamma_stat_")) {
    const RooAbsReal* nomVar= dynamic_cast<const RooAbsReal*>(ws->arg("nom_"+name));
    if (!nomVar) return 0;
    Double_t nom= nomVar->getVal();
    if (const RooAbsReal* sigVar= dynamic_cast<const RooAbsReal*>(ws->arg(name+"_sigma"))) {
      nomVal= nom;
      nomErr= sigVar->getVal();
    } else {
      if (nom<=0.0) return 0;
      nomVal= 1.0;
      nomErr= sqrt(1.0/nom);
      nomErr += 0.5/nom; //Poisson has slightly longer tails than Gaussian. This correction approximately matches p-value at +5sigma
    }
    return 2;
  } else {
    RooRealVar* glob= ws->var(name+"_In");
    if (!glob)  glob= ws->var("nom_"+name);
    if (!glob)  glob= ws->var("nom_alpha_"+name);
    if (!glob)  glob= ws->var("RNDM_"+name);
    if (!glob) return 0;
    // Can't check error, so assume constrained unit Gaussian (don't use for unconstrained)
    if (glob->getVal()!=0.0) return 0;
    nomVal= 0.0;
    nomErr= 1.0;
    return 1;
  }
}

void WorkspaceCalculator::AdjustParameterRanges (RooWorkspace* ws, const RooArgSet* parms,
                                                 const RooAbsPdf* pdf, const RooArgSet* globalObservables,
                                                 double stdDevRange, double widerRange,
                                                 const char* snapshotName, int verbose, int adjustRanges,
                                                 double fixStatError)
{
  RooArgSet* unconstrained= SelectConstrained (parms, pdf, globalObservables, false);
  if (verbose>=1 && unconstrained->getSize()>0) {
    cout << "Unconstrained parameters:-"<<endl;
    unconstrained->Print("v");
  }
  RooArgSet* origParms= 0;
  if (snapshotName && *snapshotName) {
    origParms= dynamic_cast<RooArgSet*>(parms->snapshot());
    ws->loadSnapshot (snapshotName);
  }
  RooArgSet setConst;
  for (RooLinkedListIter it= parms->iterator(); TObject* o= it.Next();) {
    RooRealVar* v= dynamic_cast<RooRealVar*>(o);
    if (!v) continue;
    TString name= v->GetName();
    RooRealVar* v0= origParms ? dynamic_cast<RooRealVar*>(origParms->find(name)) : v;
    if (v0->isConstant()) continue;
    Double_t lo= v->getMin(), hi= v->getMax();
    int mod= 0;
    int type=0; // 0=unknown, 1=Gaussian constraint, 2=gamma_stat Poisson or Gaussian, 3=unconstrained
    static const char* typeName[]= {"unknown      ", "  constrained", "gamma        ", "unconstrained"};
    Double_t val= 0.0, err= -1.0, range= stdDevRange;
    if (unconstrained->find(name)) type= 3;
    else                           type= NominalValue (ws, name, val, err);
    if (type==0) {
      cout << "Can't classify NP "; v->Print();
      continue;
    }
    if (type==0 || type==3) {
      val= v->getVal();
      err= v->getError();
      range= widerRange;
    }
    if (type==2 && fixStatError!=0.0 && err>0.0 && (fixStatError>0 ? (err<fixStatError) : (err>-fixStatError))) {
      if (v->isConstant()) continue;
      setConst.addOwned(*new RooRealVar(v->GetName(),v->GetTitle(),val));
      cout << "Set "<<name<<" constant at "<<v->getVal()<< " ("<<Form("%.2f",100.0*err)<<"% error)" << endl;
      continue;
    }
    if (adjustRanges<=0) continue;
    if (range<0.0) err= range;
    Double_t nerr= (err>0.0) ? range*err : err;
    static WildcardList normParms ("ATLAS_norm_* ATLAS_Norm_* alpha_ATLAS_Norm_*");
    if (normParms.Match (name) && lo>0.0) {
      mod++;
      v->setMin(0.0);
    } else if (err>0.0 && lo < (val-nerr)) {
      mod++;
      v->setMin(val-nerr);
    }
    static WildcardList expandParms ("slope_* gamma_stat_*");
    if (err>0.0 && ((adjustRanges>=2 && expandParms.Match(name)) ||
                    hi > (val+nerr))) {
      mod++;
      v->setMax(val+nerr);
    }
    if (mod && verbose>=0) {
      cout << "Adjust "<<typeName[type]<<" NP "<<name<<" range from ("<<lo<<" - "<<hi<<") to ("<<v->getMin()<<" - "<<v->getMax()<<")"<<endl;
    } else if (verbose>=1) {
      cout << "Keep   "<<typeName[type]<<" NP "<<name<<" range "<<lo<<" - "<<hi<<endl;
    }
  }
  if (origParms) {
    *const_cast<RooArgSet*>(parms)= *origParms;  // doesn't restore original limits
    delete origParms;
  }
  delete unconstrained;
  for (RooLinkedListIter it= setConst.iterator(); RooRealVar* vc= dynamic_cast<RooRealVar*>(it.Next());) {
    RooRealVar* v= dynamic_cast<RooRealVar*>(parms->find(vc->GetName()));
//  v->setVal(vc->getVal());
    v->setConstant();
  }
}


int WorkspaceCalculator::SetupWorkspace()
{
   Info("SetupWorkspace","Running HypoTestCalculator #%d with test statistic #%d on the workspace '%s'",
        calculatorType,testStatType,ws->GetName());

   // Load model and data from workspace

   RooMsgService::instance().setGlobalKillBelow (RooFit::MsgLevel(max (int(RooFit::DEBUG),
                                                                  min (int(RooFit::FATAL),
                                                                  RooFit::PROGRESS-verbose))));

   if (optimize & kLoadInitialSnapshot) LoadSnapshot (ws, initialSnapshot, "Loaded initial snapshot '%s'");

   if (verbose>=3) ws->Print();

   // get the modelConfig out of the file
   sbModel = dynamic_cast<RooStats::ModelConfig*>(ws->genobj(modelSBName));

   if (!sbModel) {
      Error("SetupWorkspace","ModelConfig '%s' does not exist",modelSBName.Data());
      return 0;
   }
   // check the model
   if (!sbModel->GetPdf()) {
      Error("SetupWorkspace","Model '%s' has no PDF",modelSBName.Data());
      return 0;
   }
   if (!sbModel->GetParametersOfInterest()) {
      Error("SetupWorkspace","Model '%s' has no POI",modelSBName.Data());
      return 0;
   }
   if (!sbModel->GetObservables()) {
      Error("SetupWorkspace","Model '%s' has no observables ",modelSBName.Data());
      return 0;
   }

   allParams= new RooArgSet (*sbModel->GetParametersOfInterest(), *sbModel->GetNuisanceParameters());

   if (!wsEdit.empty()) {
     Info("SetupWorkspace","Edit workspace '%s':",ws->GetName());
     TString oldPdfName= sbModel->GetPdf()->GetName();
     TString newPdfName= oldPdfName+"_edit";
     for (size_t i=0, n=wsEdit.size(); i<n; i++) {
       TString e= wsEdit[i].c_str();
       bool hadOld= e.Contains("OLDPDF");
       e.ReplaceAll ("OLDPDF",oldPdfName);
       e.ReplaceAll ("NEWPDF",newPdfName);
       if (verbose>=0) cout << e << endl;
       RooAbsArg* a= ws->factory(e);
       if (verbose>=0 && a) { cout << "-> "; a->Print(); }
       if (hadOld) {
         RooAbsPdf* newPdf= dynamic_cast<RooAbsPdf*>(a);
         if (newPdf && newPdf->GetName() != oldPdfName) {
           sbModel->SetPdf(newPdf->GetName());
           cout << "Set " << newPdf->GetName() << " as PDF in model " << sbModel->GetName() << endl;
         }
       }
     }
   }

   if (optimize & kFixCache) {
     RooArgSet* comps= sbModel->GetPdf()->getComponents();
     for (RooFIter iter= comps->fwdIterator(); RooAbsArg* a=iter.next();) {
       if (a->InheritsFrom("RooStarMomentMorph")) {
         int error= 0;
//       if (RooStarMomentMorph* p= dynamic_cast<RooStarMomentMorph*>(a)) p->fixCache();
         a->Execute("fixCache","",&error);  // ((RooStarMomentMorph*)a)->fixCache() but doesn't need RooStarMomentMorph compiled in
         if (error) cout << "Error "<<error<<" fixing cache for " << a->ClassName() << "::"<< a->GetName() << endl;
         else if (verbose>=1) cout << "Fixed cache for " << a->ClassName() << "::"<< a->GetName() << endl;
       }
     }
     delete comps;
   }

   if (optimize & kInterpolateErrors) ModifyInterpolationForAll (ws, interpCodeNorm, interpCodeShape);

   RooArgList setvars;
   if (dataName=="ASIMOV" || (dataName.Index("=")!=kNPOS && !*Scan (dataName, setvars) && setvars.getSize()>0)) {
     RooArgSet* snap= SetParameters (&setvars, false);
     RooArgSet* snapvars= dynamic_cast<RooArgSet*>(ws->components().selectCommon(*snap));
     if (verbose>=1 && snapvars->getSize()>0) {
       cout << "Generate Asimov dataset with:" << endl;
       snapvars->Print("v");
     } else if (verbose>=0)
       cout << "Generate Asimov dataset" << endl;
     RooStats::AsymptoticCalculator::SetPrintLevel (verbose>=1 ? verbose : verbose+1);  // -q->0 def,-v->1, -vv->2

     RooArgSet asimovGlobs;
     data= RooStats::AsymptoticCalculator::MakeAsimovData (*sbModel, RooArgSet(), asimovGlobs);
     //     data= RooStats::AsymptoticCalculator::GenerateAsimovData (*sbModel->GetPdf(), *sbModel->GetObservables());
     *snapvars = *snap;
     delete snapvars;
     delete snap;

     if (verbose>=0) {
       cout << "Created Asimov dataset ";
       if (verbose>=1) {
         cout << Form("with %d/%g entries: ", data->numEntries(), data->sumEntries());
         data->Print("v");
       } else
         data->Print();
     }

     resetGlobsSnapshot= "asimovData_initialGlobs";
     ws->saveSnapshot("asimovData_initialGlobs",asimovGlobs);
     const char* asimovGlobsName= "asimovData_globObs";
     cout << Form("Save Asimov global observables snapshot '%s' in workspace '%s'", asimovGlobsName, ws->GetName()) << endl;
     ws->saveSnapshot (asimovGlobsName, asimovGlobs, kTRUE);
     ws->loadSnapshot (asimovGlobsName);

     const char* asimovName= "asimovData";
     if (ws->data(asimovName))  // don't overwrite existing asimovData dataset if we write out this workspace
       kept.Add(data);
     else {
       cout << Form("Save Asimov dataset '%s' in workspace '%s'", asimovName, ws->GetName()) << endl;
       data->SetName(asimovName);
       ws->import(*data);
       delete data;
       data= ws->data(asimovName);
     }
   } else {
     data = ws->data(dataName);
     if (!data) {
       Error("SetupWorkspace","Dataset '%s' does not exist",dataName.Data());
       return 0;
     } else if (verbose>=0) {
       cout << "Using dataset ";
       data->Print();
     }
     if (optimize & kVectorStore) data->convertToVectorStore();
   }

   if (noSystematics) {
     // case of no systematics
     // remove nuisance parameters from model
     const RooArgSet* nuisPar = sbModel->GetNuisanceParameters();
     if (nuisPar && nuisPar->getSize() > 0) {
       Info("SetupWorkspace","Switch off all systematics by setting them constant to their initial values");
       const_cast<RooArgSet*>(nuisPar)->setAttribAll("Constant",kTRUE);
     }
   }

   if (optimize & kSetConst) {
     if (const RooArgSet* nuisPar = sbModel->GetNuisanceParameters()) {
       for (RooLinkedListIter it = nuisPar->iterator(); TObject* o = it.Next();) {
         RooRealVar* p= dynamic_cast<RooRealVar*>(o);
         if (!p) continue;
         if (TString(p->GetName()).Contains("gamma_stat_")) continue;
         //           if (TString(p->GetName()).BeginsWith("ATLAS_EM_")) continue;
         p->setConstant();
         if (verbose>=0) {
           cout << "Set nuisance parameter constant: "; p->Print();
         }
       }
     }
   }

   if (optimize & (kAdjustRanges|kFixStatError)) {
     AdjustParameterRanges (ws, sbModel->GetNuisanceParameters(), sbModel->GetPdf(), sbModel->GetGlobalObservables(),
                            newParameterRanges, 3.0*newParameterRanges, "conditionalNuis_muhat", verbose,
                            !(optimize&kAdjustRanges) ? 0 : (optimize&kExpandRanges) ? 2 : 1,
                            (optimize&kFixStatError) ? fixStatErrorMin : 0.0);
   }

   if (optimize & kGenerateBinned) {
     if (RooDataSet* newData= SetDatasetBinning (sbModel->GetPdf(), data, 500, "BINNEDGEN",
                                                 "PttEtaConvVBFCat* combCat_PttEtaConvVBFCat* \
                                                  *_gg_201? *_gg_1112 *_gg_201?_es *_gghc_201? *_gghc_201?_es \
                                                  *_gges_ht_201[12] *_gges_ht_1112 \
                                                  unconv_* conv_* 7TeV_unconv_* 7TeV_conv_* MoriondCat?_201? MoriondCat??_201? *_yy *_yy_ATLAS *_gg *_gg_atlas",
                                                 "ATLAS_H_* combCat_ATLAS_H_*",
                                                 verbose)) {
       if (optimize & kSaveWorkspace) {
         dataName= newData->GetName();
         if (ws->import (*newData))
           cerr << "Error importing dataset "<<newData->GetName()<<endl;
         else {
           RooAbsData* d= ws->data(dataName);
           if (d) {
             delete newData;
             data= d;
           } else {
             cerr << Form("Newly imported dataset %s could not be retrieved",dataName.Data()) << endl;
             data= newData;
             calcObjs.Add(data);
           }
         }
       } else {
         data= newData;
         calcObjs.Add(data);
       }
     }
   }

   if (optimize & kBinnedLikelihoodOpt) {
     // Activate binned likelihood calculation for binned models
     static WildcardList H4l ("ATLAS_Signal_*");  // Don't use for unbinned H->4l components
     for (RooFIter iter= ws->components().fwdIterator(); RooAbsArg* arg= iter.next();) {
       if (arg->IsA() == RooRealSumPdf::Class() && !H4l.Match(arg->GetName()) && !arg->getAttribute("BinnedLikelihood")) {
         arg->setAttribute ("BinnedLikelihood");
         if (verbose >= 1) {
           cout << "Enable BinnedLikelihood for ";
           arg->Print();
         }
       }
     }
   }

   if ((optimize & kSetGlobalObservable) &&
       sbModel->GetGlobalObservables() && sbModel->GetGlobalObservables()->getSize() > 0 &&
       !sbModel->GetPdf()->getStringAttribute("DefaultGlobalObservablesTag")) {
     if (verbose >= 0) cout << "Set GLOBAL_OBSERVABLE attribute on global observables" << endl;
     sbModel->GetPdf()->setStringAttribute("DefaultGlobalObservablesTag","GLOBAL_OBSERVABLE");
     const_cast<RooArgSet*>(sbModel->GetGlobalObservables())->setAttribAll("GLOBAL_OBSERVABLE");
   }

   if (optimize & kAddPseudoData)
     AddPseudoData (data, sbModel->GetPdf(),
                    "ch_*_201[12] *Cat_*201[12] combCat_*Cat_*201[12]_llllcat_1112",
                    "100:180:0.1,180:600:1",
                    verbose);

   if (optimize & kDisableNumInt) {
     static const WildcardList niPdfs ("ATLAS_H_*"); // HCP WS doesn't need this: "modelunc_ATLAS_H_* *_llll_201?"
     RooArgSet pdfs= ws->allPdfs();
     for (RooLinkedListIter it = pdfs.iterator(); TObject* o= it.Next();) {
       if (!niPdfs.Match(o->GetName())) continue;
       RooRealSumPdf* sumPdf= dynamic_cast<RooRealSumPdf*>(o);
       if (!sumPdf) continue;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,6)
       if (!sumPdf->getForceNumInt()) continue;  // nothing to do, so skip message
#endif
       if (verbose>=0) cout << "Disable numeric integration for PDF "<<sumPdf->GetName()<<endl;
       sumPdf->forceNumInt(false);
     }
   }

   SetParameters (varSettings);

   if (optimize & kAddErrors) AddErrors (sbModel, verbose>=1);

   return 1;
}


RooArgSet* WorkspaceCalculator::SetParameters (const RooAbsCollection* setvars, bool updateNP)
{
   // Updates sbModel->GetNuisanceParameters() if updateNP (default), otherwise returns snapshot of original parameter values.
   const RooArgSet* poi=   sbModel->GetParametersOfInterest();
   RooArgSet        np=   *sbModel->GetNuisanceParameters();
   RooArgSet        vars= *sbModel->GetParametersOfInterest();
   vars.add              (*sbModel->GetNuisanceParameters());
   RooArgSet        gobs= *sbModel->GetGlobalObservables();

   RooArgSet *allVars=0, *oldSnap=0, *constraints=0, *constVars=0;
   int nnp=0, ngl=0;
   for (RooLinkedListIter it= setvars->iterator(); const RooRealVar* v= dynamic_cast<const RooRealVar*>(it.Next());) {
     // select exact name in workspace, or by wildcard in POI+NPs, or else by wildcard in workspace.
     // can also specify %ATTRIB to select by attribute.
     RooAbsArg* fv= dynamic_cast<RooAbsArg*>(ws->arg(v->GetName()));
     if (!dynamic_cast<RooRealVar*>(fv) && !dynamic_cast<RooCategory*>(fv)) fv=0;
     RooAbsCollection* sel= 0;
     if (fv)
       sel= new RooArgSet(*fv);
     else if (*v->GetName() == '%') {
       sel= vars.selectByAttrib(v->GetName()+1,kTRUE);
       if (!sel->getSize()) {
         if (!allVars) {
           allVars= new RooArgSet(ws->allVars());
           allVars->add(ws->allCats());
         }
         delete sel;
         sel= allVars->selectByAttrib(v->GetName()+1,kTRUE);
       }
     } else {
       sel= vars.selectByName(v->GetName());
       if (!sel->getSize()) {
         if (!allVars) {
           allVars= new RooArgSet(ws->allVars());
           allVars->add(ws->allCats());
         }
         delete sel;
         sel= allVars->selectByName(v->GetName());
       }
     }
     if (const char* excl= v->getStringAttribute("excludeNames")) {
       RooAbsCollection* rem= (*excl == '%') ? sel->selectByAttrib(excl+1,kTRUE) : sel->selectByName(excl);
       sel->remove(*rem);
       delete rem;
     }
     if (!updateNP) {
       RooArgSet* snap= dynamic_cast<RooArgSet*>(sel->snapshot());
       if (!oldSnap) oldSnap= snap;
       else          oldSnap->addOwned(*snap,kFALSE);
     }
     int nsel= 0;
     for (RooLinkedListIter jt= sel->iterator(); TObject* no= jt.Next();) {
       RooAbsArg* nv= dynamic_cast<RooAbsArg*>(no);
       RooRealVar  *nr=          dynamic_cast<RooRealVar* >(nv);
       RooCategory *nc= nr ? 0 : dynamic_cast<RooCategory*>(nv);
       if (!nr && !nc) continue;
       nsel++;
       int add=0, modc=0;
       RooAbsCollection* osnap= RooArgSet(*nv).snapshot();
       if (nv->isConstant() != v->isConstant()) {
         modc++;
         if        (v->isConstant()) {
           if (np.remove(*nv,kTRUE,kTRUE)) add=-1;
         } else if (!poi->find(*nv)) {
           if (np.add   (*nv,kTRUE))       add= 1;
         }
         if (add) nnp++;
         if (nr) nr->setConstant (v->isConstant());
         else    nc->setConstant (v->isConstant());
       } else if (!nv->isConstant() && !vars.find(*nv)) {
         if (np.add (*nv,kTRUE)) {
           add= 1;
           nnp++;
         }
       }
       int modg=0, modv=0, modl=0, mode=0;
       RooArgSet* modgl=0;
       if (updateNP && add && nr && (optimize & kUpdGlobalObservables)) {
         if (!constraints) {  // do slow getAllConstraints only once
           RooArgSet* allPdfVars= sbModel->GetPdf()->getParameters(*sbModel->GetObservables());
           constVars= dynamic_cast<RooArgSet*>(allPdfVars->selectByAttrib("Constant",kTRUE));
           constraints= getAllConstraints (*sbModel->GetPdf(), *sbModel->GetObservables(), allPdfVars, verbose);
           delete allPdfVars;
         }
         RooArgSet nvs= *nv;
         modgl= findGlobalObservable (sbModel->GetPdf(), sbModel->GetObservables(), &nvs, constVars, constraints, verbose);
         if (modgl && modgl->getSize()>0) {
           if        (add==-1) {
             if (gobs.remove(*modgl))    modg=add;
           } else if (add== 1) {
             if (gobs.add(*modgl,kTRUE)) modg=add;
           }
         }
         if (modg) ngl++;
       }
       if (nr) {
         if      (!isnan(v->getMin()) && !isnan(v->getMax())) {
           if (v->getMin()!=nr->getMin() || v->getMax()!=nr->getMax()) modl++;
           nr->setRange(v->getMin(),v->getMax());
         } else if (!isnan(v->getMin())) {
           if (v->getMin()!=nr->getMin()) modl++;
           nr->setMin(v->getMin());
         } else if (!isnan(v->getMax())) {
           if (v->getMax()!=nr->getMax()) modl++;
           nr->setMax(v->getMax());
         }
         if (!isnan(v->getVal())) {
           if (v->getVal()!=nr->getVal()) modv++;
           nr->setVal(v->getVal());
         }
         if (v->hasError()) {
           if (v->getError()!=nr->getError()) mode++;
           nr->setError(v->getError());
         }
       } else {
         Double_t vv= v->getVal();
         if (!isnan(vv)) {
           Int_t vi= int(vv+0.5);
           if (fabs(vv-vi)>=1e-10 || !nc->isValidIndex(vi))
             Error("SetupWorkspace","Variable '%s' bad category index %g",nc->GetName(),vv);
           else {
             if (vi!=nc->getIndex()) modv++;
             nc->setIndex(vi);
           }
         }
         if (!isnan(v->getMin()) || !isnan(v->getMax()))
           Error("SetupWorkspace","Variable '%s' category cannot set range %g:%g",nc->GetName(),v->getMin(),v->getMax());
       }
       if (verbose>=0) {
         if      (!updateNP) cout << "Set ";
         else if (add== 1)   cout << "New NP    ";
         else if (add==-1)   cout << "Remove NP ";
         else                cout << "Set       ";
         if (nr) {
           nv->Print("i");
           if (modv || modl || modc || mode)
             if (RooRealVar* ov= dynamic_cast<RooRealVar*>(osnap->find(*nv))) {
               cout << "(was ";
               if      (modv) cout << ov->getVal() << (ov->isConstant() ? "C" : "");
               else if (modc) cout << (ov->isConstant() ? "const" : "free");
               if ((modv||modc) && (modl|mode)) cout << ' ';
               if (mode) {
                 if (ov->hasError()) cout << "+/- " << ov->getError();
                 else                cout << "+/- NONE";
               }
               if (mode && modl) cout << ' ';
               if (modl) cout << "L(" << (ov->hasMin()?Form("%g",ov->getMin()):"-INF") << " - " << (ov->hasMax()?Form("%g",ov->getMax()):"+INF") << ")";
               cout << ") ";
             }
           if (modg) cout << "and " << (modgl->getSize()>1 ? Form("%d ",modgl->getSize()) : "") << "GObs " << modgl->contentsString() << (modgl->getSize()>1 ? " ***" : "");
         } else {
           cout << nc->IsA()->GetName() << "::" << nc->GetName() << " = " << nc->getLabel() << "(idx = " << nc->getIndex() << ") " << (nc->isConstant() ? "C " : "");  // skip extra \n
           if (modv || modc)
             if (RooCategory* ov= dynamic_cast<RooCategory*>(osnap->find(*nv))) {
               cout << "(was ";
               if      (modv) cout << ov->getIndex() << (ov->isConstant() ? "C" : "");
               else if (modc) cout << (ov->isConstant() ? "const" : "free");
               cout << ")";
             }
         }
         cout << endl;
       }
       delete modgl;
       delete osnap;
     }
     delete sel;
     if (!nsel) Error("SetupWorkspace","Variable '%s' does not exist in workspace",v->GetName());
   }
   if (updateNP) {
     if (nnp>0) sbModel->SetNuisanceParameters(np);
     if (ngl>0) sbModel->SetGlobalObservables(gobs);
   }
   delete allVars;
   delete constraints;
   delete constVars;
   if (!updateNP) {
     if (!oldSnap) oldSnap= new RooArgSet;
     oldSnap->setName(Form("%s_snapshot",setvars->GetName()));
   }
   return oldSnap;
}



int WorkspaceCalculator::SetupMinimizer()
{
   // Set up minimizer

   if ((optimize & kMinuit2) && minimizerType.IsNull()) minimizerType = "Minuit2";
   ROOT::Math::MinimizerOptions::Default("Minuit2").SetValue("StorageLevel",0);  // don't save intermediate results (Hessian can use lots of RAM)
   if      (optimize & kStrategy2) ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
   else if (optimize & kStrategy0) ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   if (minimizerType.IsNull()) minimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
   else                                        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerType);
   if (fitTol>=0.0) ROOT::Math::MinimizerOptions::SetDefaultTolerance(fitTol);
   if (maxFunctionCalls>0) ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(maxFunctionCalls);
   if (fitPrec>=0.0) ROOT::Math::MinimizerOptions::SetDefaultPrecision(fitPrec*ROOT::Minuit2::MnMachinePrecision().Eps());
   ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(verbose);
   if (verbose>=1) ROOT::Math::MinimizerOptions::PrintDefault(minimizerType);
   return 1;
}

int WorkspaceCalculator::SetupModel()
{
   if (poiName.size()>0) {
     RooArgSet  newPOIs;
     RooArgList oldPOIs= *sbModel->GetParametersOfInterest();
     int modpoi= labs(int(poiName.size())-int(oldPOIs.getSize()));
     for (size_t i=0, n=poiName.size(); i<n; i++) {
       RooAbsArg* v= ws->arg(poiName[i].c_str());
       if (!v) {
         Error("SetupModel","No POI '%s'",poiName[i].c_str());
         return 0;
       }
       newPOIs.add(*v);
       RooAbsArg* vo= oldPOIs.at(i);
       if (vo && v->namePtr()!=vo->namePtr()) modpoi++;
     }
     if (modpoi) {
       int modnp= 0;
       RooArgSet np= *sbModel->GetNuisanceParameters();
       RooArgSet rempoi= *sbModel->GetParametersOfInterest();
       rempoi.remove (newPOIs);
       if (rempoi.getSize()>0) {
         if (verbose>=0) {
           cout << "Remove POIs:" << endl;
           rempoi.Print("v");
         }
         RooArgSet* rempoiFree= dynamic_cast<RooArgSet*>(rempoi.selectByAttrib("Constant",kFALSE));
         np.add (*rempoiFree,kTRUE);
         if (rempoiFree->getSize()>0) modnp++;
         delete rempoiFree;
       }
       RooArgSet addpoi= newPOIs;
       addpoi.remove (*sbModel->GetParametersOfInterest());
       if (addpoi.getSize()>0) {
         if (verbose>=0) {
           cout << "Add POIs:" << endl;
           addpoi.Print("v");
         }
         np.remove (addpoi,kTRUE);
         if (addpoi.getSize()>0) modnp++;
       }
       if (modnp) sbModel->SetNuisanceParameters(np);
       sbModel->SetParametersOfInterest(newPOIs);
     }
   } else {
     for (RooLinkedListIter it= sbModel->GetParametersOfInterest()->iterator(); TObject* o= it.Next();)
       poiName.push_back (o->GetName());
   }

   if (sbModel->GetParametersOfInterest()->getSize() <= 0) {
     Error("SetupModel","No POIs defined");
     return 0;
   }

   // Create snapshots of the POI for S+B and B-only models
   if (sbPOI.empty()) {
     if (sbModel->GetParametersOfInterest()->getSize()==1 && !isnan(sbPOIdefault))
       sbPOI.push_back(sbPOIdefault);
     else {
       for (RooLinkedListIter it= sbModel->GetParametersOfInterest()->iterator(); TObject* o= it.Next();)
         if (RooRealVar* v= dynamic_cast<RooRealVar*>(o))
           sbPOI.push_back (v->getVal());
     }
   }
   const RooArgSet* sbSnapshot = SetValueInSnapshot (sbModel, sbPOI);

   // Don't print ModelConfig, because this loads the snapshot, resetting the POI
   // if (verbose>=1) sbModel->Print();

   if (modelBName.Length()>0) {
     bModel = dynamic_cast<RooStats::ModelConfig*>(ws->genobj(modelBName));
     if (!bModel || bModel == sbModel) {
       Error("SetupModel","The background model '%s' does not exist",modelBName.Data());
       return 0;
     }
     if (bModel->GetParametersOfInterest()->getSize() <= 0) {
       Error("SetupModel","No B-only POIs defined");
       return 0;
     }
     if (noSystematics) {
       if (const RooArgSet* bnuisPar = bModel->GetNuisanceParameters())
         if (bnuisPar && bnuisPar->getSize() > 0) const_cast<RooArgSet*>(bnuisPar)->setAttribAll("Constant",kTRUE);
     }
   } else {
     modelBName = Form("%sBkg",modelSBName.Data());
     Info("SetupModel","Copy background model '%s' from '%s'",modelBName.Data(),modelSBName.Data());
     bModel = dynamic_cast<RooStats::ModelConfig*>(sbModel->Clone(modelBName));
     modelBName= bModel->GetName();
     calcObjs.Add(bModel);
   }

   if (bPOI.empty()) {
     if (bModel->GetParametersOfInterest()->getSize()==1 && !isnan(bPOIdefault))
       bPOI.push_back(bPOIdefault);
     else {
       for (RooLinkedListIter it= bModel->GetParametersOfInterest()->iterator(); TObject* o= it.Next();)
         if (RooRealVar* v= dynamic_cast<RooRealVar*>(o))
           bPOI.push_back (v->getVal());
     }
   }
   const RooArgSet* bSnapshot = SetValueInSnapshot (bModel, bPOI);

   // check model has global observables when there are nuisance pdf
   // for the hybrid case the globobs are not needed
   if (calculatorType != 1) {
      bool hasNuisParam = (sbModel->GetNuisanceParameters() && sbModel->GetNuisanceParameters()->getSize() > 0);
      bool hasGlobalObs = (sbModel->GetGlobalObservables()  && sbModel->GetGlobalObservables() ->getSize() > 0);
      if (hasNuisParam && !hasGlobalObs) {
         // try to see if model has nuisance parameters first
         RooAbsPdf* constrPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisanceConstraintPdf_sbmodel");
         if (constrPdf) {
            Warning("SetupModel","Model %s has nuisance parameters but no global observables associated - the effect of the nuisance parameters will not be treated correctly",sbModel->GetName());
         }
      }
   }

   sbPdf = sbModel->GetPdf();
    bPdf =  bModel->GetPdf();

#ifdef USE_Analytic
   if ((optimize & kAnalytic) && sbPdf == bPdf) {
     if (RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(sbPdf)) {
       Info("SetupModel","Replace PDF with HistFactorySimultaneous");
       RooStats::HistFactory::HistFactorySimultaneous* hsPdf = new RooStats::HistFactory::HistFactorySimultaneous (*simPdf);
       hsPdf->Print();
       sbPdf= hsPdf;
        bPdf= hsPdf;
       calcObjs.Add(hsPdf);
       if (optimize & kSaveWorkspace) {
         if (verbose>=0) cout << "Replace PDF "<<hsPdf->GetName()<<" in model" << endl;
         sbModel->SetPdf(*hsPdf);
          bModel->SetPdf(*hsPdf);
       }
     }
   }
#endif

   if (BkgIsAlt()) {
     nullModel    = sbModel;
      altModel    =  bModel;
     nullPdf      = sbPdf;
      altPdf      =  bPdf;
     nullSnapshot = sbSnapshot;
      altSnapshot =  bSnapshot;
     nullPOI      = sbPOI;
      altPOI      =  bPOI;
   } else {
     nullModel    =  bModel;
      altModel    = sbModel;
     nullPdf      =  bPdf;
      altPdf      = sbPdf;
     nullSnapshot =  bSnapshot;
      altSnapshot = sbSnapshot;
     nullPOI      =  bPOI;
      altPOI      = sbPOI;
   }

   if (nullMLESnapshot.Length()==0 && nullPOI.size()>0)
     nullMLESnapshot= "conditionalNuis_"+Join(nullPOI,"_");
   if ( altMLESnapshot.Length()==0 &&  altPOI.size()>0)
      altMLESnapshot= "conditionalNuis_"+Join( altPOI,"_");

   // Don't print ModelConfig, because this loads the snapshot, resetting the POI
   // if (verbose>=1) bModel->Print();
   if (verbose>=0) {
     cout << "Null/"  << nullModelName << " snapshot " << Join(*nullSnapshot,", ",1)
          << ", Alt/" <<  altModelName << " snapshot " << Join( *altSnapshot,", ",1) << endl;
   }

   return 1;
}


int WorkspaceCalculator::SetupInitialParms()
{
   // Set up the parameter of interest

   poimin.resize (sbModel->GetParametersOfInterest()->getSize(), NaN);
   poimax.resize (sbModel->GetParametersOfInterest()->getSize(), NaN);
   size_t i= 0;
   for (RooLinkedListIter it = sbModel->GetParametersOfInterest()->iterator(); TObject* o = it.Next(); i++) {
     RooRealVar* poi= dynamic_cast<RooRealVar*>(o);
     if (!poi) continue;
     if (isnan(poimax[i]))       poimax[i]= poi->getMax();
     else                        poi->setMax(poimax[i]);
     if (poimin.size()==1 && isnan(poimin[i])) {
       if      (testStatType==5) poimin[i]= OneSidedPositiveMin;
       else if (testStatType==6) poimin[i]= -poimax[i];
     }
     if (isnan(poimin[i]))       poimin[i]= poi->getMin();
     else                        poi->setMin(poimin[i]);
   }
   if (verbose>=0) {
     cout << "POI range";
     if (poimin.size()==1) {
       cout << ' ';
       sbModel->GetParametersOfInterest()->first()->Print();
     } else {
       cout << ":-" << endl;
       sbModel->GetParametersOfInterest()->Print("v");
     }
   }

   bool doFit = initialFit;
   // Need initial fit for LEP test stat or adaptive importance sampling, but maybe snapshot is OK
   if (initialFit == -1 && !(optimize & kInitialFit)) doFit = false;
   //   if (initialFit == -1 && (testStatType != 0 && samplerType != 0 && !(optimize & kInitialFit))) doFit = false;

   if (doFit) DoInitialFit();

   // print a message in case of LEP test statistics because it affects result by doing or not doing a fit
   if (testStatType == 0) {
      if (!doFit)
         Info("SetupInitialParms","Using LEP test statistic - an initial fit is not done and the TS will use the nuisances at the model value");
      else
         Info("SetupInitialParms","Using LEP test statistic - an initial fit has been done and the TS will use the nuisances at the best fit value");
   }

   // If we didn't fit, assume snapshot has best fit
   if (RooRealVar* poi = dynamic_cast<RooRealVar*>(sbModel->GetParametersOfInterest()->first())) {
     poihat  = poi->getVal();
     poierr  = poi->getError();
   }

   // Update allParams, so it includes all parameters that were or now are in the POI/NP lists.
   // We want to be sure that it can reset any vars that might be in a snapshot we might load.
   allParams->add (*sbModel->GetParametersOfInterest());
   allParams->add ( *bModel->GetParametersOfInterest());
   allParams->add (*sbModel->GetNuisanceParameters());
   allParams->add ( *bModel->GetNuisanceParameters());

   ws->saveSnapshot("initialNuisAndPoi",*allParams);

   if (optimize & kFastInit) {
     // Faster init, without splitting NPs into constrained, unconstrained, and gamam_stat.
     nullParams= new RooArgSet;
     nullParams->add (*sbModel->GetParametersOfInterest());
     nullParams->add (*sbModel->GetNuisanceParameters());
     minosParams= new RooArgSet (*nullParams);
   } else {
     RooArgSet* allVars= nullPdf->getParameters(*data);
     RooArgSet* vars= dynamic_cast<RooArgSet*>(allVars->selectByAttrib("Constant",kFALSE));
     delete allVars;
     vars->remove(*nullModel->GetParametersOfInterest());
     vars->sort();
     nullParams= new RooArgSet (*nullModel->GetParametersOfInterest());
     nullParams->add (*vars);

     minosParams= new RooArgSet (*nullModel->GetParametersOfInterest());
     if (!constrainedOnly) minosParams->add (*vars);

     RooArgSet* constrained= SelectConstrained (vars, nullModel->GetPdf(), nullModel->GetGlobalObservables());
     // H->bb alpha_Systtbar3JNorm is constrained, but doesn't look it (see mail from Gabriel Facini, 19/11/2013)
     RooAbsCollection* really_constrained= vars->selectByName("alpha_Systtbar3JNorm,alpha_Systtbar3JNorm_*,ATLAS_Hbb_Systtbar3JNorm_*");
     constrained->add(*really_constrained,kTRUE);
     delete really_constrained;
     vars->remove(*constrained);
     RooAbsCollection* gamma= vars->selectByName("*gamma_stat_*");
     vars->remove(*gamma);

     if (constrainedOnly) minosParams->add (*constrained);

     if (verbose >= 0) {
       cout << "====== POIs ======" << endl;
       nullModel->GetParametersOfInterest()->Print("v");
       if (constrained->getSize() > 0) {
         cout << "====== Constrained NPs ======"              << endl;
         constrained->Print("v");
       }
       if (vars->getSize() > 0) {
         cout << "====== Unconstrained NPs ======"            << endl;
         vars->Print("v");
       }
       if (gamma->getSize() > 0) {
         cout << "====== Histogram bin statistics NPs ======" << endl;
         gamma->Print("v");
       }
     }
     delete gamma;
     delete constrained;
     delete vars;
   }

   if ((optimize & (kMinosData|kAltMinos)) || !minosSetNames.empty() || minosSetSize>0) {
     if (!minosSetNames.empty()) {
       minosSet= new RooArgSet;
       for (size_t j= 0, n=minosSetNames.size(); j<n; j++) {
         RooAbsArg* v= nullParams->find(minosSetNames[j].c_str());
         if (v) minosSet->add(*v);
         else   cout << "MINOS parameter " << minosSetNames[j] << " not found" << endl;
       }
     } else
       minosSet= new RooArgSet (*minosParams);
     int nset= 0;
     if (minosSetSize>0) {
       RooArgList vlist = *minosSet;
       int n= vlist.getSize();
       int iset1= minosSetNum*minosSetSize, iset2= iset1+minosSetSize;
       nset= (n+minosSetSize-1)/minosSetSize;
       if (iset2>n) iset2= n;
       delete minosSet;
       minosSet= new RooArgSet;
       for (int iset= iset1; iset<iset2; iset++) minosSet->add (vlist[iset]);
     }
     if (minosSet->getSize()<=0) {
       cout << "Empty MINOS parameter set - skip MINOS" << endl;
       delete minosSet; minosSet= 0;
     }
     if (verbose>=0) {
       if (minosSet && minosSet->getSize() < nullParams->getSize()) {
         if (minosSetSize>0) cout << "Run MINOS for parameter set "<<minosSetNum<<"/"<<nset<<":";
         else                cout << "Run MINOS for parameters";
         for (RooFIter it= minosSet->fwdIterator(); const RooAbsArg* a= it.next();)
           cout << ' ' << a->GetName();
         cout << endl;
       } else if (minosSet && (optimize & (kMinosData|kAltMinos))) {
         cout << "Run MINOS for all parameters" << endl;
       }
     }
   }
   return 1;
}


void WorkspaceCalculator::DoInitialFit()
{
   // fit the data first (need to use constraint)

   // do the fit : By doing a fit the POI snapshot (for S+B) is set to the fit value
   // and the nuisance parameters nominal values will be set to the fit value.
   // This is relevant when using LEP test statistics

   RooRealVar* poi = 0;
   if (sbModel->GetParametersOfInterest()->getSize()==1) {
     poi = dynamic_cast<RooRealVar*>(sbModel->GetParametersOfInterest()->first());
     Double_t mu_start = poi->getMin() + 1.05*poi->getError();
     if (poi->getVal() < mu_start && mu_start+1.05*poi->getError() < poi->getMax()) poi->setVal(mu_start);
     Info("DoInitialFit","Doing a first fit to the observed data. Initial %s = %g",poi->GetName(),poi->getVal());
   } else {
     Info("DoInitialFit","Doing a first fit to the observed data. Initial POI:-");
     sbModel->GetParametersOfInterest()->Print("v");
   }
   RooArgSet constrainParams;
   if (sbModel->GetNuisanceParameters()) constrainParams.add(*sbModel->GetNuisanceParameters());
   RooArgSet* constrainParamsFree= dynamic_cast<RooArgSet*>(constrainParams.selectByAttrib("Constant",kFALSE));
   TStopwatch tw;
   RooFitResult* fitres = sbPdf->fitTo (*data,
//                                      RooFit::InitialHesse(false),   // default=false
                                        RooFit::Hesse(false),
                                        RooFit::Minimizer(minimizerType, "Migrad"),
                                        RooFit::Strategy(0),
                                        RooFit::Verbose(verbose>=2),
                                        RooFit::PrintLevel(verbose-2),
                                        RooFit::Constrain(*constrainParamsFree),
                                        RooFit::SumW2Error(kFALSE),
                                        RooFit::Save(true));
   if (fitres->status() != 0) {
     Warning("DoInitialFit","Fit to the model failed - try with strategy 1 and perform first an Hesse computation");
     fitres = sbPdf->fitTo (*data,
                            RooFit::InitialHesse(true),
                            RooFit::Hesse(false),
                            RooFit::Minimizer(minimizerType, "Migrad"),
//                          RooFit::Strategy(1),  // default=1
                            RooFit::Verbose(verbose>=2),
                            RooFit::PrintLevel(verbose-2),
                            RooFit::Constrain(*constrainParamsFree),
                            RooFit::SumW2Error(kFALSE),
                            RooFit::Save(true));
   }
   delete constrainParamsFree;
   if (fitres->status() != 0)
     Warning("DoInitialFit"," Fit still failed - continue anyway.....");
   if (verbose>=0) {
     cout << "Time for fitting: "; tw.Print();
   }

   if (poi)
     Info("DoInitialFit","Best Fit value : %s  = %g +/- %g",poi->GetName(), poi->getVal(), poi->getError());
}


template <class PLTS>
PLTS* WorkspaceCalculator::SetupPLTS (PLTS* profll, const char* tsname)
{
  // Common setup between ProfileLikelihoodTestStat and ProfileLikelihoodTestStatEnhanced.

  profll->SetPrintLevel(verbose+1);
#ifdef USE_DetailedOutput
  if      (detailedOutput>=1)             profll->EnableDetailedOutput(1,(detailedOutput>=2));
#ifdef USE_ToyMCImportanceSampler
  else if (doStatsTree && samplerType>=0) profll->EnableDetailedOutput(2);  // still need weights
#endif
#endif

  if (sbModel->GetConditionalObservables() && sbModel->GetConditionalObservables()->getSize() > 0) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,2)
    Info("SetupPLTS","Set %d conditional observables",sbModel->GetConditionalObservables()->getSize());
    if (verbose>=1) sbModel->GetConditionalObservables()->Print("v");
    profll->SetConditionalObservables(*sbModel->GetConditionalObservables());
#else
    Warning("SetupPLTS","Cannot define %d conditional observables in %s",sbModel->GetConditionalObservables()->getSize(),profll->GetVarName().Data());
#endif
  }
  if        (testStatType == 2) {
    Info("SetupPLTS","Using two-sided %s",tsname);
  } else if (testStatType == 3) {
    profll->SetOneSided(true);
    Info("SetupPLTS","Using one-sided %s",tsname);
#ifdef ROOT53201
  } else if (testStatType == 5) {
    profll->SetOneSidedDiscovery();
    Info("SetupPLTS","Using one-sided discovery %s",tsname);
  } else if (testStatType == 6) {
    profll->SetOneSidedDiscovery();
    profll->SetSigned();
    Info("SetupPLTS","Using signed discovery %s",tsname);
  } else if (testStatType == 7) {
    profll->SetOneSided();
    profll->SetSigned();
    Info("SetupPLTS","Using signed %s",tsname);
#endif
  } else {
    delete profll; profll= 0;
  }
  if (profll) {
    profll->SetMinimizer(minimizerType);
    if      (optimize & kStrategy2) profll->SetStrategy(2);
    else if (optimize & kStrategy0) profll->SetStrategy(0);
    if (fitTol>=0.0) profll->SetTolerance(fitTol);
  }
  return profll;
}


#ifdef USE_ProfileLikelihoodTestStatEnhanced
void WorkspaceCalculator::SetInitialSnapshotTS (RooStats::ProfileLikelihoodTestStatEnhanced* profll)
{
  TString snapName, snapNameCond;
  Int_t j= initialSnapshotTS.Index(":",0);
  if (j!=kNPOS) {
    snapName=     initialSnapshotTS (0, j);
    snapNameCond= initialSnapshotTS (j+1, 999999);
  } else
    snapName= initialSnapshotTS;
  RooArgSet* allVars= nullPdf->getParameters(*data);
  RooArgSet* vars= dynamic_cast<RooArgSet*>(allVars->selectByAttrib("Constant",kFALSE));
  delete allVars;
  RooArgSet *initSnap=0, *initSnapCond=0;
  RooArgSet setvars, setvarsCond;
  if (snapName.Index("=")!=kNPOS && !*Scan (snapName, setvars) && setvars.getSize()>0) {
    cout << "Use variables to initialise unconditional fit:" << endl;
    setvars.Print("v");
    initSnap= dynamic_cast<RooArgSet*>(setvars.snapshot());
  } else
    initSnap= GetSnapshot (ws, snapName, "Use snapshot '%s' to initialise unconditional fit", vars);
  if (snapNameCond.Length()>0) {
    if (snapNameCond.Index("=")!=kNPOS && !*Scan (snapNameCond, setvarsCond) && setvarsCond.getSize()>0) {
      cout << "Use variables to initialise conditional fit:" << endl;
      setvarsCond.Print("v");
      initSnapCond= dynamic_cast<RooArgSet*>(setvarsCond.snapshot());
    } else
      initSnapCond= GetSnapshot (ws, snapNameCond, "Use snapshot '%s' to initialise conditional fit", vars);
  }
  delete vars;
  profll->SetInitialValues (initSnap, initSnapCond);
  calcObjs.Add(initSnap);
  if (initSnapCond) calcObjs.Add(initSnapCond);
}

void WorkspaceCalculator::SetInitGlobObs (RooStats::ProfileLikelihoodTestStatEnhanced* profll)
{
  const RooArgSet *nuisPar, *globObs;
  if ((nuisPar = sbModel->GetNuisanceParameters()) && (globObs = sbModel->GetGlobalObservables())) {
    RooConstVar* one= new RooConstVar ("1", "1", 1.0);
    RooArgSet* initSet= new RooArgSet();
    calcObjs.Add(one);
    calcObjs.Add(initSet);
    for (RooLinkedListIter it = nuisPar->iterator(); TObject* o = it.Next();) {
      RooRealVar* np= dynamic_cast<RooRealVar*>(o);
      if (!np) continue;
      TString name= np->GetName();
      static WildcardList noInit ("\
        ATLAS_B_EFF ATLAS_CLUSTER ATLAS_EMB ATLAS_E_EFF ATLAS_JER_TOP ATLAS_JES ATLAS_JES_TOP \
        ATLAS_M_SCALE ATLAS_T_EFF ATLAS_VH_B4_EFF ATLAS_VH_C_EFF ATLAS_VH_JER ATLAS_VH_JES \
        ATLAS_VH_WBBPTW ATLAS_VH_WTOPNORM1 ATLAS_VH_ZBBPTZ ATLAS_VH_ZTOPNORM1 FakeRate UEPS \
        alpha_ATLAS_MU_RESCALE_lvlv_2012 alpha_ATLAS_Z_LINESHAPE_SYS_lvlv \
        alpha_ATLAS_Z_SCALEF_EF_mm0jBK_lvlv alpha_AddOn_lh alpha_CLUSTERBKLH_lh alpha_EtoTau_lh \
        alpha_GenAccq2Z_lh alpha_JESBKGO_ll alpha_JESBKLH_lh alpha_PILEUPBKLH_lh \
        alpha_QCDMethod_hh alpha_SameSign_lh alpha_TESBKHH_hh alpha_T_SF_ll alpha_ZHHMethod_hh \
        alpha_fakes_ll alpha_pdf_qqbarlh_lh alpha_pdf_qqbarll_ll lumi gamma_stat_chan0j_bin_10_ll \
        gamma_stat_chan0j_bin_9_ll gamma_stat_chan1j_bin_10_ll gamma_stat_chan2jVBF_bin_1_ll \
        gamma_stat_chan2jVH_bin_2_ll gamma_stat_channel*_lh \
      ");
      if (noInit.Match(name)) continue;
      RooAbsArg* glob= globObs->find(name+"_In");
      if (!glob) glob= globObs->find("nom_"+name);
      if (!glob) glob= globObs->find("RNDM_"+name);
      if (!glob) continue;
      RooAbsArg* sf= one;
      if (name.Contains("gamma_stat_")) {
        string npName= name.Data(), tauName= npName+"_tau";
        for (size_t upos= string::npos;;) {
          // cout << "Guessed: " << tauName << endl;
          if (RooAbsReal* tauVar= dynamic_cast<RooAbsReal*>(ws->arg(tauName.c_str()))) {
            if (tauVar->getVal()==0.0) {
              cout << "bad tau parameter: ";
              tauVar->Print();
              break;
            }
            string invName= tauName+"_inverse";
            sf= new RooConstVar (invName.c_str(), invName.c_str(), 1.0/tauVar->getVal());
            calcObjs.Add(sf);
            break;
          }
          upos = npName.rfind('_',upos);
          if (upos==0 || upos == string::npos) {
            cerr << "no tau parameter for " << name << endl;
            break;
          }
          tauName = npName.substr(0,upos) + "_tau";
          upos--;
        }
        if (sf==one) continue;
      }
      initSet->addOwned(*new RooProduct (name, name+"_initFit", RooArgList(*glob,*sf)));
    }
    if (initSet->getSize()>0) {
      Info ("SetInitGlobObs","Initialise fit using %d global observables",initSet->getSize());
      if (verbose>=2) initSet->Print("v");
      profll->SetInitialValues(initSet);
    }
  }
}
#endif


int WorkspaceCalculator::SetupTestStat()
{
   // build test statistics

   if (!(optimize & kReuseNLL)) {
     RooStats::ProfileLikelihoodTestStat::SetAlwaysReuseNLL(kFALSE);
#ifdef USE_ProfileLikelihoodTestStatEnhanced
     RooStats::ProfileLikelihoodTestStatEnhanced::SetAlwaysReuseNLL(kFALSE);
#endif
     RooStats::SimpleLikelihoodRatioTestStat::SetAlwaysReuseNLL(kFALSE);
     RooStats::RatioOfProfiledLikelihoodsTestStat::SetAlwaysReuseNLL(kFALSE);
   }

   if        (testStatType == 0) {
     RooStats::SimpleLikelihoodRatioTestStat* slrts =
       new RooStats::SimpleLikelihoodRatioTestStat (*nullPdf,*altPdf);
     slrts->SetNullParameters (RooArgSet(*nullSnapshot,*sbModel->GetNuisanceParameters()));
     slrts->SetAltParameters  (RooArgSet( *altSnapshot, *bModel->GetNuisanceParameters()));
     testStat = slrts;
   } else if (testStatType == 1 || testStatType == 11) {
     // ratio of profile likelihood - need to pass snapshot for the alt
     RooStats::RatioOfProfiledLikelihoodsTestStat* ropl =
       new RooStats::RatioOfProfiledLikelihoodsTestStat (*nullPdf, *altPdf, altSnapshot);
     ropl->SetSubtractMLE(testStatType==11);
     ropl->SetPrintLevel(verbose+1);
     ropl->SetMinimizer(minimizerType);
     if      (optimize & kStrategy2) ropl->SetStrategy(2);
     else if (optimize & kStrategy0) ropl->SetStrategy(0);
     testStat = ropl;
   } else if (testStatType == 2 || testStatType == 3 || testStatType == 5 || testStatType == 6 || testStatType == 7) {
#ifdef USE_ProfileLikelihoodTestStatEnhanced
     if (optimize & kPLTSEnhanced) {
       RooStats::ProfileLikelihoodTestStatEnhanced* profll = new RooStats::ProfileLikelihoodTestStatEnhanced (*nullPdf);
       profll->SetPrintLevel(verbose+1);
       profll->FlagData(data);  // clear in SetupCalculator, if used
       profll->SkipToys(skipToys);
       profll->PlotParameterProfiles(parmProf);
       if (fitTol0>=0.0) profll->SetTolerance0(fitTol0);
       if (optimize & kSkipDataTS)     profll->SetOptimize("D",0);
       if (optimize & kSaveInfo)       profll->SetOptimize("Sr",0);
//     if (optimize & kGenerateBinned) profll->SetOptimize("B",0);
       profll->SetOptimize(optimizeTS);
       if      (optimize & kUseSavedToys) profll->UseSavedToyData(ws);
       else if (optimize & kSaveInfo)     profll->SaveToyData (ws, runToys ? nullModel->GetGlobalObservables() : 0);
       else                               profll->SaveToyData (0,  runToys ? nullModel->GetGlobalObservables() : 0);
       if (minosSet && (optimize & kMinosData)) {
         if (minosSet->getSize() < nullParams->getSize())
           profll->UseMinos(minosSet);
         else
           profll->UseMinos();
       }

       if (doStatsTree) {
         profll->FillTree();
         if (optimize & kAddSLRTS) {
           RooStats::SimpleLikelihoodRatioTestStat* slrts = new RooStats::SimpleLikelihoodRatioTestStat (*nullPdf,*altPdf);
           slrts->SetNullParameters (*nullSnapshot);
           slrts->SetAltParameters  ( *altSnapshot);
           profll->SetComparisonTestStat(slrts);
         }
       }

       if (initialSnapshotTS.Length()>0) SetInitialSnapshotTS (profll);
       if (optimize & kInitGlobObs)      SetInitGlobObs       (profll);
       if (optimize & kFixInitially) {
         RooArgSet* constrained= SelectConstrained (nullModel->GetNuisanceParameters(), nullPdf, nullModel->GetGlobalObservables());
         if (constrained->getSize()>0) {
           profll->SetFixInitially(constrained);
           calcObjs.Add(constrained);
         } else
           delete constrained;
       }
       testStat = SetupPLTS (profll,"profile likelihood test statistic (enhanced)");
     } else
#endif
       testStat = SetupPLTS (new RooStats::ProfileLikelihoodTestStat (*nullPdf));
   } else if (testStatType == 4) {
     if (nullModel->GetParametersOfInterest()->getSize()!=1) {
       Warning ("SetupTestStat","MaxLikelihoodEstimateTestStat only works with one POI. We have:-");
       nullModel->GetParametersOfInterest()->Print("v");
     }
     RooRealVar* poi = dynamic_cast<RooRealVar*>(nullModel->GetParametersOfInterest()->first());
     RooStats::MaxLikelihoodEstimateTestStat* maxll = new RooStats::MaxLikelihoodEstimateTestStat (*nullPdf, *poi);
     testStat = maxll;
#ifdef USE_ProfiledLikelihoodRatioTestStatExt
   } else if (testStatType == 8) {
     ProfiledLikelihoodTestStatOpt* profll = new ProfiledLikelihoodTestStatOpt (*nullModel->GetObservables(), *nullPdf, nullModel->GetNuisanceParameters(), *nullSnapshot, RooArgList(), RooArgList(), verbose+1);
     testStat= profll;
#endif
   } else if (testStatType == 8) {
     testStat= new RooStats::NumEventsTestStat;
   }
   if (!testStat) {
     Error("SetupTestStat","Invalid - test statistic type = %d supported values are only : 0 (SLR) , 1 (Tevatron) , 2 (PLR), 3 (PLR1), 4(MLE)"
#if defined(USE_ProfileLikelihoodTestStatEnhanced) || defined(ROOT53201)
           ", 5 (PLR1D), 6 (PLR+-D), 7 (PLR+-)"
#endif
           ,testStatType);
     return 0;
   }
   Info("SetupTestStat","Using %s test statistic",testStat->GetVarName().Data());
   return 1;
}


int WorkspaceCalculator::SetupSampler()
{
   // Set up toy MC sampler
   if (calculatorType == 2 || calculatorType == 3) return 1;

   if (samplerType<0) toymcs = new RooStats::ToyMCSampler (*testStat, 100);
#ifdef USE_ToyMCImportanceSampler
   else {
     RooStats::ToyMCImportanceSampler* toymcis= new RooStats::ToyMCImportanceSampler (*testStat, 100);
     if (!(optimize & kReuseNLL)) toymcis->SetReuseNLL(false);
     toymcis->SetExpIncreasingNumToysPerDensity();  // use more toys for importance densities with larger mu
//   toymcis->SetEqualNumToysPerDensity();          // or not
     int numImpSnaps=0;
     RooRealVar* poi = dynamic_cast<RooRealVar*>(nullModel->GetParametersOfInterest()->first());
#ifdef USE_ToyMCImportanceSampler_profile
     if (RooArgSet* condMLE= GetSnapshot (ws, "conditionalNuis_muhat",
                                          Form("Use '%%s' snapshot to generate weighted toys at muhat = %g",poihat),
                                          sbModel->GetNuisanceParameters())) {
       if (samplerType==0 && nStdDevOverlap>0.0)
         numImpSnaps= toymcis->CreateImpDensitiesForOnePOIAdaptively (*sbPdf, *sbModel->GetParametersOfInterest(), *poi,
                                                                      nStdDevOverlap, bPOI[0], condMLE, data);
       else {
         numImpSnaps= toymcis->CreateNImpDensitiesForOnePOI          (*sbPdf, *sbModel->GetParametersOfInterest(), *poi,
                                                                      samplerType,    bPOI[0], condMLE, data);
       }
       delete condMLE;
     }
#else
     if (samplerType==0 && nStdDevOverlap>0.0)
       numImpSnaps= toymcis->CreateImpDensitiesForOnePOIAdaptively (*sbPdf, *sbModel->GetParametersOfInterest(), *poi,
                                                                    nStdDevOverlap, bPOI[0]);
     else
       numImpSnaps= toymcis->CreateNImpDensitiesForOnePOI          (*sbPdf, *sbModel->GetParametersOfInterest(), *poi,
                                                                    samplerType,    bPOI[0]);
#endif
     Info("SetupSampler","Create %d sets of weighted samples",numImpSnaps);
     toymcs= toymcis;
   }
#else
   else {
     Error("SetupSampler","Invalid - test statistic sampler type = %d",samplerType);
     return 0;
   }
#endif
   Info("SetupSampler","Using %s",toymcs->IsA()->GetName());

#ifdef USE_DetailedOutput
   if (optimize & kAddSLRTS) {
     RooStats::SimpleLikelihoodRatioTestStat* slrts = 0;
     slrts = new RooStats::SimpleLikelihoodRatioTestStat(*nullPdf,*altPdf);
     slrts->SetNullParameters (*nullSnapshot);
     slrts->SetAltParameters  ( *altSnapshot);
     if (detailedOutput) slrts->EnableDetailedOutput();
     toymcs->AddTestStatistic (slrts);
     testStat2= slrts;
   }
#endif

   if (useNumberCounting<0)
     useNumberCounting= (!sbPdf->canBeExtended() &&
                         ! bPdf->canBeExtended() &&
                         data->numEntries()==1) ? 1 : 0;
   if (useNumberCounting) {
     Info("SetupSampler","Use number counting");
     toymcs->SetNEventsPerToy(1);
   }
// if (optimize & kGenerateBinned) {   // always enabled in case BINNEDGEN already set in workspace. kGenerateBinned enables setting of the attribute.
     toymcs->SetGenerateBinnedTag("BINNEDGEN");
     toymcs->SetGenerateAutoBinned();
// }

   if (optimize & kUseMultiGen) {
     if (samplerType>=0)
       Warning("SetupSampler","MultiGen not enabled - does not work with importance sampling");
     else if (calculatorType == 1 && (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters()))
       Warning("SetupSampler","MultiGen not enabled - does not work with HybridCalculator and nuisance parameters");
     else
       toymcs->SetUseMultiGen(kTRUE);
   }

   return 1;
}


int WorkspaceCalculator::SetupCalculator()
{
   // build hypotest calculator for running the inverter

   if        (calculatorType == 0) {
     RooStats::FrequentistCalculator* fc = new RooStats::FrequentistCalculator (*data, *altModel, *nullModel, toymcs);
     fc->SetToys(ntoys,ntoysAlt);
#ifdef USE_StoreFitInfo
     if (detailedOutput) fc->StoreFitInfo();
#endif
#ifdef USE_SetConditionalMLEs
     RooArgSet *nullMLE=0, *altMLE=0;
     if (optimize & kUseCondSnapshot) {
       // speed up: take MLEs from conditionalNuis snapshots in workspace
       if (ntoys>0    &&
           (nullMLE= GetSnapshot (ws, nullMLESnapshot, "Use '%s' snapshot to generate null toys",
                                  nullModel->GetNuisanceParameters()))) {
         calcObjs.Add(nullMLE);
         fc->SetConditionalMLEsNull (nullMLE);
       }
       if (ntoysAlt>0 &&
           ( altMLE= GetSnapshot (ws,  altMLESnapshot, "Use '%s' snapshot to generate alt toys",
                                   altModel->GetNuisanceParameters()))) {
         calcObjs.Add (altMLE);
         fc->SetConditionalMLEsAlt  (altMLE);
       }
     }
     // If we aren't generating any toys for null or alt, then we can use any old snapshot
     // to skip the profiling.
     if (ntoys<=0    && !nullMLE) {
       nullMLE= dynamic_cast<RooArgSet*>(nullModel->GetNuisanceParameters()->snapshot());
       calcObjs.Add(nullMLE);
       fc->SetConditionalMLEsNull (nullMLE);
     }
     if (ntoysAlt<=0 &&  !altMLE) {
        altMLE= dynamic_cast<RooArgSet*>( altModel->GetNuisanceParameters()->snapshot());
       calcObjs.Add (altMLE);
       fc->SetConditionalMLEsAlt  (altMLE);
     }
#endif
     calc= fc;
   } else if (calculatorType == 1) {
     RooStats::HybridCalculator*      hc = new RooStats::HybridCalculator      (*data, *altModel, *nullModel, toymcs);
     hc->SetToys(ntoys,ntoysAlt);

     // remove global observables from ModelConfig (this is probably not needed anymore in 5.32)
     nullModel->SetGlobalObservables(RooArgSet());
      altModel->SetGlobalObservables(RooArgSet());

     // check for nuisance prior pdf in case of nuisance parameters
     if (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters()) {
       RooAbsPdf * nuisPdf = 0;
       if (nuisPriorName.Length()>0) nuisPdf = ws->pdf(nuisPriorName);
       // use prior defined first in bModel (then in SbModel)
       if (!nuisPdf)  {
         Info("SetupCalculator","No nuisance pdf given for the HybridCalculator - try to deduce pdf from the model");
         if (bPdf && bModel->GetObservables())
           nuisPdf = RooStats::MakeNuisancePdf(*bModel,"nuisancePdf_bmodel");
         else
           nuisPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisancePdf_sbmodel");
       }
       if (!nuisPdf) {
         Error("SetupCalculator","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified");
         return 0;
       }
       const RooArgSet * nuisParams = (bModel->GetNuisanceParameters()) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
       RooArgSet * np = nuisPdf->getObservables(*nuisParams);
       if (np->getSize() == 0) {
         Warning("SetupCalculator","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
       }
       delete np;
       hc->ForcePriorNuisanceAlt (*nuisPdf);
       hc->ForcePriorNuisanceNull(*nuisPdf);
     }
     calc= hc;
   } else if (calculatorType == 2
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,32,1)
              || calculatorType == 3) {
     RooStats::AsymptoticCalculator::SetPrintLevel (verbose);
     RooStats::AsymptoticCalculator*  ac = new RooStats::AsymptoticCalculator  (*data, *altModel, *nullModel, (calculatorType==3));
#else
             ) {
     RooStats::AsymptoticCalculator::SetPrintLevel (verbose);
     RooStats::AsymptoticCalculator*  ac = new RooStats::AsymptoticCalculator  (*data, *altModel, *nullModel);
#endif
     if      (testStatType == 3) ac->SetOneSided(true);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,4)
     else if (testStatType == 5) ac->SetOneSidedDiscovery(true);
#endif
     calc= ac;
   } else {
     Error("SetupCalculator","Invalid calculator type = %d - supported values are: 0 (Frequentist), 1 (Hybrid), 2 (Asymptotic), 3 (Asymptotic+Asimov)",
           calculatorType);
     return 0;
   }
   Info("SetupCalculator","Using %s",calc->IsA()->GetName());
   if (calculatorType != 2 && calculatorType != 3)
     Info("SetupCalculator","Generate %d %s toys and %d %s toys",ntoys,nullModelName,ntoysAlt,altModelName);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,4)
   if (optimize & kReuseAltToys) calc->UseSameAltToys();
#endif
#ifdef USE_ProfileLikelihoodTestStatEnhanced
   if (RooStats::ProfileLikelihoodTestStatEnhanced* profll= dynamic_cast<RooStats::ProfileLikelihoodTestStatEnhanced*>(testStat)) {
     profll->SetHypoTestCalculator(calc);
     profll->FlagData(0);  // no longer needed
   }
#endif
   return 1;
}

int WorkspaceCalculator::RunCalculatorProof()
{
  if (optimize & kNoFits) {
    Info("RunCalculatorProof","SKIP RUNNING FITS");
    return 1;
  }
  // can speed up using proof-lite
  RooStats::ProofConfig* pc = 0;
  if (toymcs && nworkers > 1) {
    pc = new RooStats::ProofConfig (*ws, nworkers, "", kFALSE);
    toymcs->SetProofConfig(pc);    // enable PROOF
  }
  int ok= RunCalculator();
  delete pc;
  return ok;
}

int WorkspaceCalculator::GetTestStatInfo()
{
#ifdef USE_ProfileLikelihoodTestStatEnhanced
   if (RooStats::ProfileLikelihoodTestStatEnhanced* profll= dynamic_cast<RooStats::ProfileLikelihoodTestStatEnhanced*>(testStat)) {
     TList& info= profll->GetInfo();
     static WildcardList objnames ("fitted_poi fit_result fit_result_cond minos_results");
     TObjArray selobj(objnames.NumTerms());
     bool save_fitted_poi= false;
     for (TIter it= &info; TObject* o= it();) {
       Int_t ind= objnames.MatchInd(o->GetName());
       if (ind>=0) selobj[ind]= o;   // select fitted_poi and fitres. For fitres, only use last try
       if (!(ind>=0 && (optimize & kSaveWorkspace))) {
         if (*it.GetOption()) {
           saved.Add (o, it.GetOption());
           if (strcmp (it.GetOption(), "fitted_poicov") == 0) save_fitted_poi= true;
         } else
           saved.Add (o);
       }
     }
     // keep fitted_poi if have we fitted_poicov, so we know which parameter is which (fitted_poicov is a TMatrixDSym).
     if (save_fitted_poi && !saved.FindObject("fitted_poi"))
       if (TObject* o= info.FindObject("fitted_poi")) saved.Add (o);

     RooAbsCollection* fitted_poi= dynamic_cast<RooAbsCollection*>(selobj[0]);
     if (fitted_poi)
       if (RooRealVar* firstPOI= dynamic_cast<RooRealVar*>(fitted_poi->first())) {
         poihat= firstPOI->getVal();
         poierr= firstPOI->getError();
       }

     if (optimize & kSaveWorkspace) {
       const char*       ucmles= "ucmles";
       const char* uncondNPsnap= "conditionalNuis_muhat";
       if (const RooFitResult* ures= dynamic_cast<const RooFitResult*>(selobj[1])) {
         RooArgSet uvars= ures->floatParsFinal();
         if (verbose>=0) cout << Form("Save %d values from RooFitResult %s to snapshot '%s'",
                                      uvars.getSize(), ures->GetName(), ucmles) << endl;
         if (verbose>=2) uvars.Print("v");
         ws->saveSnapshot (ucmles, uvars, true);
         if (fitted_poi) {
           uvars.remove (*fitted_poi, false, true);
           if (verbose>=0) cout << Form("Save %d values from RooFitResult %s to snapshot '%s'",
                                        uvars.getSize(), ures->GetName(), uncondNPsnap) << endl;
           if (verbose>=2) uvars.Print("v");
           ws->saveSnapshot (uncondNPsnap, uvars, true);
         }
       }
       if (const RooFitResult* cres= dynamic_cast<const RooFitResult*>(selobj[2])) {
         RooArgSet cvars= cres->floatParsFinal();
         if (verbose>=0) cout << Form("Save %d values from RooFitResult %s to snapshot '%s'",
                                      cvars.getSize(), cres->GetName(), nullMLESnapshot.Data()) << endl;
         if (verbose>=2) cvars.Print("v");
         ws->saveSnapshot (nullMLESnapshot, cvars, true);
       }
       if (RooArgSet* minos= dynamic_cast<RooArgSet*>(selobj[3])) {
         if (verbose>=0) cout << Form("Save %d values from RooArgSet %s to snapshot '%s'",
                                      minos->getSize(), minos->GetName(), minos->GetName()) << endl;
         if (verbose>=2) minos->Print("v");
         ws->saveSnapshot (minos->GetName(), *minos, true);
         minos_results= minos;
       }
       if (!saved.FindObject(ws)) saved.Add(ws);
     }

     // don't let info.Clear() delete objects that are now in saved
     for (TIter it= &saved; TObject* o= it();) info.Remove(o);
     info.Clear();

     statsTree= profll->GetTree(true);
     if (statsTree) saved.Add (statsTree);
     profll->FillTree(false);
   }
#endif

#ifdef USE_StoreFitInfo
   if (calc)
     if (const RooArgSet* fitInfo= calc->GetFitInfo()) {
       RooAbsCollection* fitInfoSnap = fitInfo->snapshot();
       fitInfoSnap->setName("fitInfo");
       saved.Add(fitInfoSnap);
     }
#endif
   return 1;
}


#ifdef USE_ToyMCImportanceSampler
Int_t
WorkspaceCalculator::GetWeights (RooAbsData* ds, const RooStats::SamplingDistribution* samp, Int_t is_alt,
                                 std::vector<Double_t>& tsval, std::vector<Double_t>& wgt,
                                 std::vector<Int_t>& alt, std::vector<Int_t>& label)
{
  if (!ds || ds->numEntries()<=0 || !samp) return 0;
  TString varname;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
  varname.Form("%s_TS0",samp->GetName());
#else
  varname= samp->GetVarName();
#endif
  const RooRealVar*  ntentry= dynamic_cast<RooRealVar*> (ds->get()->find (Form("%s_ntentry",varname.Data())));
  if (!ntentry) return 0;
  const RooRealVar*  dsvar=   dynamic_cast<RooRealVar*> (ds->get()->find (varname));
  const RooCategory* labvar=  dynamic_cast<RooCategory*>(ds->get()->find ("densityLabel"));
  Int_t nset=0, nerr=0;
  for (Int_t i=0, n=ds->numEntries(), nent= wgt.size(); i<n; i++) {
    ds->get(i);
    Int_t ent= int(ntentry->getVal()+0.5);
    if (i>=0 && i<nent && alt[ent]==-1) {
      if (dsvar)  tsval[ent]= dsvar->getVal();
      if (labvar) label[ent]= labvar->getIndex();
      alt[ent]= is_alt;
      wgt[ent]= ds->weight();
      nset++;
    } else {
      if (!nerr) cerr << "Bad "<<ntentry->GetName()<<" numbers:";
      cerr << ' ' <<ntentry->getVal();
      nerr++;
    }
  }
  if (nerr) cerr << endl;
  return nset;
}


void WorkspaceCalculator::AddWeights (TTree*& tree, RooStats::HypoTestResult* r)
{
  if (!tree) return;
  Long64_t nent= tree->GetEntries();
  if (nent<=0) return;
  if (!r) return;
  vector<Double_t> tsval(nent,NaN), wgt(nent);
  vector<Int_t> label(nent,-2), alt(nent,-1);
  Int_t nset= 0;
  nset += GetWeights (r->GetNullDetailedOutput(), r->GetNullDistribution(), 0, tsval, wgt, alt, label);
  nset += GetWeights (r->GetAltDetailedOutput(),  r->GetAltDistribution(),  1, tsval, wgt, alt, label);
  if (!detailedOutput && !(optimize & kAddSLRTS)) {  // remove output if we only needed to get weights
    delete r->GetNullDetailedOutput();
    delete r->GetAltDetailedOutput();
    r->SetNullDetailedOutput(0);
    r->SetAltDetailedOutput(0);
    r->SetAllTestStatisticsData(0);
  }
  if (!nset) {
    cerr << "Could not get weights from detailed output datasets" << endl;
    return;
  }

  Double_t pll=NaN, weight=0.0;
  Int_t sample=-2, is_alt=-1;
  tree->SetBranchAddress ("pll", &pll);
  TTree* newtree= tree->CloneTree(0);
  newtree->SetDirectory(0);
  newtree->SetAutoFlush(0);
  newtree->SetBranchAddress ("weight", &weight);
  newtree->SetBranchAddress ("is_alt", &is_alt);
  newtree->SetBranchAddress ("sample", &sample);
  for (Long64_t i=0; i<nent; i++) {
    tree->GetEntry(i);
    if (!isnan(tsval[i]) && tsval[i]!=pll) {
      cerr << "Ntuple entry "<<i<<" PLL "<<pll<<" does not match dataset value "<<tsval[i]<<endl;
    }
    weight= wgt[i];
    is_alt= alt[i];
    sample= label[i];
    newtree->Fill();
  }
  newtree->ResetBranchAddresses();
  tree->ResetBranchAddresses();
  delete tree;
  tree= newtree;
}
#endif


int WorkspaceCalculator::ReadResult (const std::vector<std::string>& fileNames)
{
  // Read results from files. Can specify multiple files with wildcards.

  TTree* limitTree= 0;
  std::map<UInt_t,string> seeds;
  bool useChain= optimize & kChainNtuple;  // since ReadTree overwrites optimize
  TFile* firstFile= 0;
  int nres= 0;
  for (size_t i= 0; i<fileNames.size(); i++) {
#ifdef USE_GLOBITER
    for (FileGlobIter files(fileNames[i].c_str()); const char* fname = files();)
#else
    for (const char* fname= fileNames[i].c_str(); fname; fname=0)
#endif
    {
      TFile* f= TFile::Open(fname);
      if (!f) {
        cerr << "Error: could not open file " << fname << endl;
        if (force) continue;
        return 0;
      }

      bool allowEmpty= force, keepFile= false;
      if (!wsInfo) {
        f->GetObject ("combWS", wsInfo);
        if (wsInfo) {
          if (verbose>=0) cout << "Read workspace '"<<wsInfo->GetName()<<"' from file " << fname << endl;
          keepFile= allowEmpty= true;
          kept.Add(f);
        }
      }

      UInt_t thisSeed=0;
      if (ReadTree (f, thisSeed) && !allowEmpty) {
        delete f;
        return 0;
      }
      if (thisSeed!=0) {
        std::pair<std::map<UInt_t,string>::iterator,bool> found= seeds.insert(std::pair<UInt_t,string>(thisSeed,fname));
        if (!found.second) {
          cerr << "WARNING: file "<<fname<<" seed "<<thisSeed<<" already used in\n         file "<<found.first->second;
          if (!force) {
            cerr << " - skip" << endl;
            delete f;
            continue;
          }
          cerr << endl;
        }
      }

      int nr= 0;
      TString resname= altResultName ? wsName : resultName;
      vector<string> resnames;
      Scan (resname, resnames);
      for (size_t j=0, n=resnames.size(); j<n; j++) {
#ifdef USE_GLOBITER
        for (TDirectoryGlobIter keys(f,resnames[j].c_str(),GlobIter::kIgnore); const char* objname = keys();)
#else
          const char* objname = resnames[j].c_str();
#endif
        {
          if (verbose>=0) cout << "Read '"<<objname<<"' from file " << fname << endl;
          int err= ReadResult (f, objname);
          if (!err) {
            nr++;
            nres++;
          } else if (err>=2) {   // could not merge
            if (force) nr++;
            else {
              delete f;
              return 0;
            }
          }
        }
      }

      if (!nr && !optionalResult) {
        cerr << "File " << fname << " does not contain a '"<<resname<<"' object" << endl;
        if (!allowEmpty) {
          delete f;
          return 0;
        }
      }

      if (!useChain) {
        TTree* sadd = 0;
        f->GetObject ("pll", sadd);
        if (sadd) {
          if (verbose>=1) cout << "File " << fname << " contains statistics ntuple '" << sadd->GetName() << "' with " << sadd->GetEntries() << " entries" << endl;
          if (doStatsTree==-1) {
            sadd->SetBranchStatus ("*", 0);
            const char* const keepvars[]= {
              "weight", "ml", "ml_cond", "pll", "is_data", "pll_init", "is_alt", "sample", "mu_gen",
              "status", "status_cond", "status_init", "tries", "tries_cond",
              "tries_retry", "fit_time", "fit_time_cond", "fit_time_retry", "toynum", "seed", "mu_fit",
              "mu_fit_init", "ml_fit", "ml_fit_cond", "ml_fit_init", "mu_init", "mu_init_err", "mu", "mu_err",
              "pllcmp", "have_uncond", "have_cond", "covqual", "covqual_cond", "covqual_init",
              "edm", "edm_cond", "edm_init"
            };
            for (size_t j=0, n=sizeof(keepvars)/sizeof(keepvars[0]); j<n; j++)
              sadd->SetBranchStatus (keepvars[j], 1);
            if (!poiName.empty())
              for (size_t j=0, n=poiName.size(); j<n; j++) {
                sadd->SetBranchStatus (poiName[j].c_str(), 1);
                sadd->SetBranchStatus (Form("%s_err",poiName[j].c_str()), 1);
              }
          }
          if (statsTree) {
            if (firstFile) {  // need to make a copy from the first file
              if (TTree* newTree = statsTree->CloneTree(0)) {
                newTree->SetDirectory(0);
                newTree->SetAutoFlush(0);
                newTree->CopyEntries(statsTree);
                delete statsTree; statsTree= newTree;
                delete firstFile; firstFile= 0;
              }
            }
            TList tadd;
            tadd.Add(sadd);
            statsTree->Merge(&tadd);
          } else if (doStatsTree!=-1 && !keepFile) {  // if this is the only file, then we don't need to read it all in
            statsTree= sadd;
            firstFile= f;
            keepFile= true;
            sadd= 0;  // don't delete it
          } else if ((statsTree = sadd->CloneTree(0))) {
            statsTree->SetDirectory(0);
            statsTree->SetAutoFlush(0);
            statsTree->CopyEntries(sadd);
          }
          if (statsTree) statsTree->ResetBranchAddresses();
          delete sadd;
          nres++;
        }
      } else if (TKey* key= f->FindKey("pll")) {
        if (verbose>=1) cout << "File " << fname << " contains statistics ntuple '" << key->GetName() << "'" << endl;
        TChain* chain;
        if (!statsTree) statsTree = chain = new TChain (key->GetName(), key->GetTitle());
        else                        chain = dynamic_cast<TChain*>(statsTree);
        chain->AddFile (fname);  // does not (re)open the file here
        nres++;
      }

      if (doStatsTree) {
        TTree* sadd = 0;
        f->GetObject ("limit", sadd);
        if (sadd) {
          if (verbose>=1) cout << "File " << fname << " contains summary ntuple '" << sadd->GetName() << "' with " << sadd->GetEntries() << " entries" << endl;
          if (limitTree) {
            TList tadd;
            tadd.Add(sadd);
            limitTree->Merge(&tadd);
          } else if ((limitTree = sadd->CloneTree(0))) {
            limitTree->SetDirectory(0);
            limitTree->SetAutoFlush(0);
            limitTree->CopyEntries(sadd);
            saved.Add(limitTree);
          }
          if (limitTree) limitTree->ResetBranchAddresses();
          delete sadd;
          nres++;
        }
      }

      RooAbsCollection* minos = 0;
      f->GetObject ("minos_results", minos);
      if (minos && minos->getSize() > 0) {
        if (verbose>=1) cout << "File " << fname << " contains " << minos->ClassName() << "::" << minos->GetName() << " with " << minos->getSize() << " entries" << endl;
        if (!minos_results) {
          minos_results= new RooArgList (minos->GetName());
          saved.Add(minos_results);
        }
        minos_results->addClone(*minos);
        delete minos;
        nres++;
      }

      if (!keepFile) delete f;
    }
  }
  if (seed==0 && seeds.size()==1) seed= seeds.begin()->first;

  if (statsTree) {
    if (doStatsTree) saved.Add(statsTree);
    else             kept .Add(statsTree);
  }
  if (firstFile) kept.Add(firstFile);
  return (nres ? 1 : 0);
}


void
WorkspaceCalculator::PrintSamplingDistribution (RooStats::SamplingDistribution* s, const char* name, int nprint)
{
  if (!s) return;
  if (s->GetSize()==0) return;
  if (s->GetSize()>=2) s->InverseCDF(0.0);  // Gives warning with weighted dist - could use instead: s->Integral(-1e20,-1e19,kFALSE); // do this first to sort GetSamplingDistribution() array
  const vector<Double_t>& vec= s->GetSamplingDistribution();
  const vector<Double_t>& wgt= s->GetSampleWeights();
  int n= vec.size();
  int nneg= 0;
  if (testStatType!=6) {
    for (int j=0; j<n; j++)
      if (vec[j]<negVal) nneg++;
  } else {
    for (int j=0; j<n; j++)
      if (vec[j]==0.0) nneg++;
  }
  cout << n << ' ' << name << " samples:";
  if (n>nprint) {
    // Print nprint samples, smallest (0) first and largest (n-1) last.
    for (int j=0; j<nprint; j++) {
      int k= int(double((n-1)*j)/(nprint-1)+0.5);
      cout << ' ' << vec[k];
      if (wgt[k]!=1.0) cout << '*' << wgt[k];
    }
    cout << "...";
  } else {
    for (int j=0; j<n;      j++) {
      cout << ' ' << vec[j];
      if (wgt[j]!=1.0) cout << '*' << wgt[j];
    }
  }
  if (nneg>0) {
    if (testStatType!=6) cout << " (" << nneg << " negative)";
    else                 cout << " (" << nneg << " zero)";
  }
  cout << endl;
}


RooStats::SamplingDistribution*
WorkspaceCalculator::ResampleToy (RooStats::SamplingDistribution* s)
{
  if (!s) return 0;
  const vector<Double_t>& vec= s->GetSamplingDistribution();
  int n= vec.size();
  if (n<=1) return 0;
  vector<Double_t> newvec(n);
  for (int i=0; i<n; i++) {
    newvec[i]= vec[RooRandom::integer(n)];
  }
  RooStats::SamplingDistribution* t=
    new RooStats::SamplingDistribution (s->GetName(), s->GetTitle(), newvec, s->GetVarName());
  delete s;
  return t;
}


RooStats::SamplingDistribution*
WorkspaceCalculator::LimitSamples (RooStats::SamplingDistribution* s, int ntoys, int skip)
{
  if (!s) return 0;
  const vector<Double_t>& vec= s->GetSamplingDistribution();
  if (vec.empty()) return 0;
  int n= vec.size();
  if (ntoys<0) ntoys= n;
  int nskip= skip>0 ? skip*ntoys : 0;
  int nnew= n-nskip;
  if      (nnew<0)     nnew= 0;
  else if (nnew>ntoys) nnew= ntoys;
  vector<Double_t> oldvec=vec, newvec(nnew);
  for (int i=0; i<n; i++) {
    int j= RooRandom::integer(oldvec.size()); // draw the same number of random numbers independently of skip and ntoys
    int k= i-nskip;
    if (k>=0 && k<nnew) newvec[k]= oldvec[j];
    oldvec.erase(oldvec.begin()+j);
  }
  if (ntoys<0) return 0;
  if (skip<=0 && ntoys>=n) return 0;
  RooStats::SamplingDistribution* t=
    new RooStats::SamplingDistribution (s->GetName(), s->GetTitle(), newvec, s->GetVarName());
  delete s;
  return t;
}


RooStats::SamplingDistribution*
WorkspaceCalculator::DropBadSamples (RooStats::SamplingDistribution* s)
{
  if (!s) return 0;
  const vector<Double_t>& vec= s->GetSamplingDistribution();
  if (vec.empty()) return 0;
  int n= vec.size();
  const vector<Double_t>& wgt= s->GetSampleWeights();
  vector<Double_t> newvec(n), newwgt(n);
  int nok=0, nnan=0, nneg=0, nlarge=0;
  for (int j=0; j<n; j++) {
    if (dropNaN && isnan(vec[j]))           { nnan++; continue; }
    if (dropNegative && vec[j]<negVal)      { nneg++; continue; }
    if (dropLarge>0.0 && vec[j]>=dropLarge) { nlarge++; continue; }
    if (wgt[j]!=1.0) haveWeights++;
    newvec[nok]= vec[j];
    newwgt[nok]= wgt[j];
    nok++;
  }
  if (nok==n) return 0;
  if (nnan) cout << "Drop "<<nnan<<" NaN samples"<<endl;
  if (nneg) cout << "Drop "<<nneg<<" negative samples"<<endl;
  if (nlarge) cout << "Drop "<<nlarge<<" samples >= "<<dropLarge<<endl;
  newvec.resize(nok);
  newwgt.resize(nok);
  RooStats::SamplingDistribution* t=
    new RooStats::SamplingDistribution (s->GetName(), s->GetTitle(), newvec, newwgt, s->GetVarName());
  delete s;
  return t;
}


RooStats::SamplingDistribution*
WorkspaceCalculator::ApplyCuts (RooStats::SamplingDistribution* s, bool use_alt)
{
  if (!s) return 0;
  if (!statsTree) return 0;
  TString hypo= (use_alt ? "Alt" : "Null");
  const Long64_t nent= statsTree->GetEntries();
  if (!nent) return 0;
  int nsamp= s->GetSize();
  if (nsamp==0) return 0;
  Int_t is_data=-1, is_alt=-1;
  Double_t pll=NaN, weight=-1.0, poi=NaN, ml=NaN, ml_cond=NaN, poicut=NaN;
  Long64_t nsel=0, ntree=0;
  TTreeFormula* select=0;
  if (cut.Length()>0) {
    select = new TTreeFormula ("selection", cut, statsTree);
    if (!select->GetNdim()) {
      cerr << "bad cut: " << cut << endl;
      delete select; select=0;
    }
  }

  statsTree->SetBranchStatus("*",0);
  SetBranchAddress (statsTree, "is_data", &is_data);
  SetBranchAddress (statsTree, "is_alt",  &is_alt);
  SetBranchAddress (statsTree, "weight",  &weight);
  bool redoPll= (overrideTestStatType>=0 &&
                 (testStatType==2 || testStatType==3 || testStatType>=5) &&
                 poiName.size()==1 && (testStatType!=5 || bPOI.size()==1) &&
                 statsTree->GetBranch("ml") &&
                 statsTree->GetBranch("ml_cond") &&
                 (testStatType==2 || statsTree->GetBranch(poiName[0].c_str())));
  if (redoPll) {
    if (testStatType!=2)
      SetBranchAddress (statsTree, poiName[0].c_str(), &poi);
    SetBranchAddress (statsTree, "ml",      &ml);
    SetBranchAddress (statsTree, "ml_cond", &ml_cond);
    if (bPOI.size()==1) poicut= bPOI[0];
  } else {
    if (!cut.Length()) {
      statsTree->ResetBranchAddresses();
      statsTree->SetBranchStatus("*",1);
      delete select;
      return 0;
    }
    SetBranchAddress (statsTree, "pll",     &pll);
  }

  if (select)
    for (Int_t i=0, n=select->GetNcodes(); i<n; i++)
      select->GetLeaf(i)->GetBranch()->SetStatus(1);
  vector<Double_t> newvec(nent), newwgt(nent);
  for (Long64_t entry=0; entry<nent; entry++) {
    pll= poi= ml= ml_cond= NaN;
    statsTree->GetEntry(entry);
    if (is_data)                         continue;
    if ((is_alt!=0) != use_alt)          continue;
    if (redoPll) {
      if (testStatType!=2 && isnan(poi)) continue;
      pll= ml_cond - ml;
      if      (testStatType==3 && poi>=poicut) pll= 0.0;
      else if (testStatType==5 && poi<=poicut) pll= 0.0;
      else if (testStatType==6 || testStatType==7) {
        if (pll<0.0) pll= 0.0;  // bad fit
        if (testStatType==6 && poi<poicut) pll= -pll;
        if (testStatType==7 && poi>poicut) pll= -pll;
      }
    }
    if (isnan(pll))                      continue;
    ntree++;
    if (select) {
      Int_t i=0, ndata=select->GetNdata();
      for (; i<ndata; i++) {
        if (select->EvalInstance(i) != 0.0) break;
      }
      if (i>=ndata) continue;
    }
    newvec[nsel]= pll;
    newwgt[nsel]= weight;
    nsel++;
  }
  statsTree->ResetBranchAddresses();
  statsTree->SetBranchStatus("*",1);
  delete select;
  if (ntree != nsamp) {
    cerr << hypo << " sampling distribution has " << nsamp
         << " entries but ntuple has " << ntree << " - don't apply cuts" << endl;
    return 0;
  }
  if (nsel == nsamp && !redoPll) return 0;
  if (verbose>=0) {
    if (nsel == nsamp)
      cout << "Convert " << ntree << " " << hypo << " entries from test statistic #"
           << overrideTestStatType << " to #" << testStatType << " (b=" << poicut << ")" << endl;
    else {
      cout << ntree-nsel << " / " << ntree << " " << hypo << " entries fail cut: "
           << cut << " - update SamplingDistribution";
      if (redoPll)
        cout << " and convert from test statistic #" << overrideTestStatType << " to #" << testStatType << " (b=" << poicut << ")";
      cout << endl;
    }
  }
  newvec.resize(nsel);
  newwgt.resize(nsel);
  RooStats::SamplingDistribution* t=
    new RooStats::SamplingDistribution (s->GetName(), s->GetTitle(), newvec, newwgt, s->GetVarName());
  delete s;
  return t;
}


int
WorkspaceCalculator::NumFailures (const RooStats::SamplingDistribution* s)
{
  if (!s) return 0;
  const vector<Double_t>& vec= s->GetSamplingDistribution();
  int n= vec.size(), nfail= 0;
  if (testStatType!=6) {
    for (int j=0; j<n; j++)
      if (vec[j]<negVal) nfail++;
  } else {
    for (int j=0; j<n; j++)
      if (vec[j]==0.0) nfail++;
  }
  return nfail;
}

Double_t WorkspaceCalculator::FillSamples (TH1* h, const RooStats::SamplingDistribution* s, bool shiftZero)
{
  if (!s) return 0;
  const vector<Double_t>& vec= s->GetSamplingDistribution();
  const vector<Double_t>& wgt= s->GetSampleWeights();
  int n= vec.size();
  Double_t nw=0.0;
  for (int j=0; j<n; j++) {
    Double_t q= 2.0*vec[j], w= wgt[j];
    if (shiftZero && q==0.0) q=-1e-14;
    if (std::isinf(q)) q= -2.0;  // HypoTestResult does not count inf in the right tail, so don't plot it there either
    h->Fill(q,w);
    nw += w;
  }
  return nw;
}

TH1* WorkspaceCalculator::NewShaded (const TH1* h, Double_t minShaded, Style_t fillStyle, Color_t fillColor)
{
  TH1* shaded= (TH1*)h->Clone(Form("%s_shaded",h->GetName()));
  shaded->SetFillColor(shaded->GetLineColor());
  shaded->SetLineWidth(0);
  if (fillStyle>=0) shaded->SetFillStyle(fillStyle);
  if (fillColor>=0) shaded->SetFillColor(fillColor);
  shaded->SetAxisRange(minShaded+0.5*shaded->GetBinWidth(1),shaded->GetXaxis()->GetXmax());
  return shaded;
}

TH1D* WorkspaceCalculator::NewTSHist (const char* name, const char* title, double ref, bool twoSided, int plotType, double negVal)
{
  double xhi=14.0;
  double xlo=twoSided?-xhi:-1.0;
  int nbin= twoSided?35:70;
  if      (plotType==2) nbin *= 10;
  else if (plotType==1) nbin= int(0.4*nbin+0.5);

  if (!isnan(ref) && 1.2*ref>xhi) xhi= min (50.0, 1.2*ref);
  double xw=xhi/nbin;
  int bin0=int(-xlo/xw);
  nbin += bin0;
  xlo= -xw*bin0;    // adjust for rounding
  if (twoSided) {   // put q0=0 in at bin centre
    xlo += 0.5*xw;
    xhi += 0.5*xw;
  } else {          // include -0.025:0 in 0 bin
    xlo += 2.0*negVal;
    xhi += 2.0*negVal;
  }
  return new TH1D (name, title, nbin, xlo, xhi);
}

void WorkspaceCalculator::wstitle()
{
  TText text;
  text.SetTextSize(0.02);
  text.SetTextFont(82);
  text.DrawTextNDC(0.01,0.01,wsfile.c_str());
}

void WorkspaceCalculator::PrintPlot (TPad* canvas)
{
  if (showCaptions) wstitle();
  if (skipPlot-- > 0) return;
  const Int_t oldErrorIgnoreLevel= gErrorIgnoreLevel;
  gErrorIgnoreLevel= kWarning;   // remove tedious info message
  canvas->Print (plotFile);
  gErrorIgnoreLevel= oldErrorIgnoreLevel;
#if defined(__CINT__) || defined(__CLING__) || defined(__ACLIC__)
  // Interactive ROOT: pause until user double-clicks or presses a key in the window before going to next plot.
  canvas->Update();
  canvas->WaitPrimitive();
#endif
  canvas->Clear();
  plotObjs.Clear();
}


Long64_t WorkspaceCalculator::Project (TTree* nt, TH1* h, const char* varexp, const char* selection,
                                       Option_t* option, Long64_t nentries, Long64_t firstentry)
{
  // Fix really stupid ROOT bug, to project onto a histogram even if TH1::AddDirectory=0.
  TDirectory *oldcwd= gDirectory, *olddir= h->GetDirectory();
  if (!olddir) h->SetDirectory(gROOT);
  h->GetDirectory()->cd();
  Long64_t n= nt->Project (h->GetName(), varexp, selection, option, nentries, firstentry);
  oldcwd->cd();
  if (!olddir) h->SetDirectory(0);
  return n;
}


TPaveText* WorkspaceCalculator::AtlasLabel (TPad* canvas, Double_t xoff, Double_t yoff, bool transparent)
{
  if (xoff>=0.0) xoff += 0.01+canvas->GetLeftMargin();
  else           xoff += 0.99-canvas->GetRightMargin();
  if (yoff> 0.0) yoff += 0.01+canvas->GetBottomMargin();
  else {
    if (yoff==0.0 && gStyle->GetPadTickY()) yoff= -0.02;
    yoff += 0.99-canvas->GetTopMargin();
  }
  Float_t siz= 0.75*gStyle->GetLabelSize();
  const int nlines= 3;
  TPaveText* label= new TPaveText (xoff, yoff-nlines*(siz+0.01), xoff+0.05, yoff, "tNDC");
  label->SetTextFont(42);
  label->SetTextSize(siz);
  label->SetTextAlign(12);
  label->SetBorderSize(0);
  label->SetFillColor(0);
  if (transparent) label->SetFillStyle(0);
  label->AddText("#scale[1.3]{#bf{#it{ATLAS}} Internal}");
  label->AddText("#sqrt{s} = 13 TeV #scale[0.6]{#int}Ldt = 3.2 fb^{-1}");
  label->Draw();
  return label;
}


TString WorkspaceCalculator::tsName (int tsType)
{
  if (tsType<0) tsType= testStatType;
  if (bPOI.size()==1) {
    if (bPOI[0]==0.0) {
      if (tsType==5 || (tsType==2 && poimin[0]==0.0)) return "q_{0}";
      if (tsType==6)                                  return "q_{0}#kern[-1]{'}";
    }
    if (tsType==3) {
      if (poimin[0]==0.0) return "#tilde{q_{#mu}}";
      else                return        "q_{#mu}";
    }
  }
  return "-2ln#Lambda";
}

TString WorkspaceCalculator::chi2Name (int tsType, int ndf)
{
  if (ndf<0) ndf=bPOI.size();
  if      (ndf==1 && asimovSig>0.0) return Form("f(#tilde{t}_{#mu}|%g)",bPOI[0]);
  else if (tsType==2) return                                    Form("#chi^{2}_{%d}", ndf);
  else                return Form("#scale[0.5]{#frac{1}{2}}#kern[0.2]{#chi^{2}_{%d}}",ndf);
}

Double_t WorkspaceCalculator::AsymFunc (const Double_t* xx, const Double_t* p)
{
  Double_t x=xx[0], N=p[0], ndf=p[1], slo=p[2], shi=p[3];
  Double_t slo2=slo*slo, shi2=shi*shi;
  Double_t f;
  if        (x <= slo2)   {
    f=       ROOT::Math::chisquared_pdf(x,ndf);
  } else if (x <= shi2) {
    f=   0.5*ROOT::Math::chisquared_pdf(x,ndf)        + ROOT::Math::gaussian_pdf(x+slo2,2.0*slo);
  } else {
    f=       ROOT::Math::gaussian_pdf(x+slo2,2.0*slo) + ROOT::Math::gaussian_pdf(x+shi2,2.0*shi);
  }
  return N*f;
}

Double_t WorkspaceCalculator::AsymPValue (const Double_t* xx, const Double_t* p)
{
  Double_t x=xx[0], N=p[0], ndf=p[1], slo=p[2], shi=p[3];
  Double_t slo2=slo*slo, shi2=shi*shi;
  Double_t f;
  if        (x <= slo2)   {
    f=       ROOT::Math::chisquared_cdf_c(x,ndf);
  } else if (x <= shi2) {
    f=   0.5*ROOT::Math::chisquared_cdf_c(x,ndf)        + ROOT::Math::gaussian_cdf_c(x+slo2,2.0*slo);
  } else {
    f=       ROOT::Math::gaussian_cdf_c(x+slo2,2.0*slo) + ROOT::Math::gaussian_cdf_c(x+shi2,2.0*shi);
  }
  return N*f;
}

Double_t WorkspaceCalculator::AsymQuantile (const Double_t* xx, const Double_t* p)
{
  // Inverse of AsymPValue calculated by binary search between p[4] and p[5]
  const Double_t eps= 1e-12;  // #iterations = -ln(eps)/ln(2) = 40
  const Double_t y=xx[0];
  Double_t xt[1];
  Double_t xlo=p[4], xhi=p[5], x= 0.5*(xhi+xlo);
  const Double_t xeps=eps*(xhi-xlo);
  while (xhi-xlo > xeps) {
    xt[0]= x;
    Double_t yt= AsymPValue (xt, p);
    if      (yt<y) xhi= x;
    else if (yt>y) xlo= x;
    else           return x;
    x= 0.5*(xhi+xlo);
  }
  return x;
}


void WorkspaceCalculator::PlotDataSet (TPad* canvas, const RooAbsPdf& pdf, const RooAbsData& dataset, bool drawComponent)
{
  // Plot dataset histograms. Based on Haoshuang Ji's asimovUtils::makePlots.

  const RooSimultaneous* sim= dynamic_cast<const RooSimultaneous*>(&pdf);
  const RooAbsCategoryLValue& cat= dynamic_cast<const RooAbsCategoryLValue&>(sim->indexCat());
  TList* datasets= dataset.split(cat, true);

  for (TIter ids= datasets; const RooAbsData* ds= dynamic_cast<const RooAbsData*>(ids());) {
    const RooAbsPdf* pdfi= sim->getPdf(ds->GetName());
    if (!pdfi) continue;  // in case dataset has extra categories, not in pdf
    RooArgSet* obs= pdfi->getObservables(ds);
    for (RooFIter ito= obs->fwdIterator(); RooAbsArg* xa= ito.next();) {
      RooRealVar* x= dynamic_cast<RooRealVar*>(xa);
      if (!x) continue;
      TString tit= Form("%s %s",ds->GetName(),x->GetName());
      Int_t nbins= x->numBins();
      if (nbins==0) nbins= 100;
      x->SetTitle(tit);

      RooPlot* plot= x->frame (RooFit::Title(tit), RooFit::Bins(nbins));
      Double_t integral= pdfi->expectedEvents(*x);

      const RooProdPdf* prodPdf;
      if (drawComponent && (prodPdf= dynamic_cast<const RooProdPdf*>(pdfi))) {
        for (RooFIter itp= prodPdf->pdfList().fwdIterator(); const RooAbsArg* va= itp.next();) {
          const RooRealSumPdf* v= dynamic_cast<const RooRealSumPdf*>(va);
          if (!v) continue;
          const RooArgList& funcList= v->funcList();
          const RooArgList& coefList= v->coefList();

          THStack* stack=  new THStack("components", "");
          TLegend* legend= new TLegend(0.63, 0.63, 0.85, 0.93);
          legend->SetFillStyle(4000);
          legend->SetBorderSize(0);
          legend->SetShadowColor(0);
          legend->SetTextFont(42);
          legend->SetTextSize(0.03);
          for (Int_t i= funcList.getSize()-1; i>=0; i--) {
            const RooProduct* func=     dynamic_cast<const RooProduct*>(funcList.at(i));
            const RooRealVar* binWidth= dynamic_cast<const RooRealVar*>(coefList.at(i));
            TH1* h= func->createHistogram (func->GetName(), *x);
            h->Scale(binWidth->getVal());
            h->SetLineColor(i+1);
            h->SetLineWidth(2);
            stack->Add(h);

            TString compName= func->GetName();
            compName.ReplaceAll("L_x_","");
            compName.ReplaceAll("ATLAS_","");
            compName.ReplaceAll("_overallSyst_x_StatUncert","");
            compName.ReplaceAll("_overallSyst_x_HistSyst","");
            compName.ReplaceAll("_overallSyst_x_Exp","");
            compName= compName(0, compName.Index("_"));
            legend->AddEntry(h, compName.Data(), "L");
          }
          plot->addObject (stack,  "same");
          plot->addObject (legend, "same");
        }
      }
      ds  ->plotOn (plot, RooFit::DataError(RooAbsData::None), RooFit::MarkerSize(1));
      pdfi->plotOn (plot, RooFit::Normalization(integral,RooAbsReal::NumEvent));

      plot->Draw();
      PrintPlot (canvas);
      delete plot;
    }
    delete obs;
    delete ds;
  }
  delete datasets;
}


TGraph* WorkspaceCalculator::PValueGraph (TF1* f)
{
  Double_t xlo,xhi;
  f->GetRange(xlo,xhi);
  Int_t n= f->GetNpx();
  Double_t xw= (xhi-xlo)/n;
  TVectorD xf(n), yf(n);
  Double_t fi=0.0;
  for (Int_t i=n-1; i>=0; i--) {
    Double_t x= xlo + xw*i;
    fi += f->Integral(x,x+xw);
    xf[i]= x;
    yf[i]= fi;
  }
  TGraph* g= new TGraph (xf,yf);
  f->TAttLine::Copy(*g);
  return g;
}

TF1* WorkspaceCalculator::AsymTF1 (int tsType, int typ, const TH1* h, Double_t norm)
{
  int nbin=100;
  Double_t xlo=0.0, xhi=25.0;
  if (h) {
    nbin=h->GetNbinsX();
    xlo=h->GetXaxis()->GetXmin();
    xhi=h->GetXaxis()->GetXmax();
  }
  Double_t xw=(xhi-xlo)/nbin;
  Int_t ndf= bPOI.size();

  TF1* f;
  if (ndf==1 && asimovSig>0.0) {
    Double_t dlo= fabs(bPOI[0]-poimin[0]);
    Double_t dhi= fabs(bPOI[0]-poimax[0]);
    Double_t slo=min(dlo,dhi)/asimovSig;
    Double_t shi=max(dlo,dhi)/asimovSig;
    if (verbose>=1) cout << Form("AsymFunc slo^2=%g, shi^2=%g",slo*slo,shi*shi)<<endl;
    if      (typ==0) f= new TF1 ("asym", AsymFunc,     0.1*xw, xhi, 6);
    else if (typ==1) f= new TF1 ("asym", AsymPValue,   0.1*xw, xhi, 6);
    else             f= new TF1 ("asym", AsymQuantile, 0.1*xw, xhi, 6);
    f->SetParameters (norm, ndf, slo, shi, xlo, xhi);
  } else {
    const char* fun=0;
    Double_t zeroBin= h ? h->GetXaxis()->GetBinLowEdge (h->GetXaxis()->FindBin (0.0)) : 0.0;
    Double_t funlo=xlo, p1=zeroBin, p3=0.0;
    if        (typ==0) {
      p1= zeroBin+xw;
      p3= 0.5 * (1.0 + (ROOT::Math::chisquared_cdf_c(zeroBin,ndf) - ROOT::Math::chisquared_cdf_c(p1,ndf)));
      if      (tsType==2) { fun=                   "[0]*ROOT::Math::chisquared_pdf(    x, [2])"; funlo= 0.1*xw; }
      else if (tsType==6) { fun=               "0.5*[0]*ROOT::Math::chisquared_pdf(abs(x),[2])"; }
      else                { fun= "x<[1] ? [3] : 0.5*[0]*ROOT::Math::chisquared_pdf(    x, [2])"; funlo= zeroBin;}
    } else if (typ==1) {
      if      (tsType==2) { fun=                   "[0]*ROOT::Math::chisquared_cdf_c(  x, [2])";  }
      else if (tsType==6) { fun=  "0.5*[0]*(1.0-sign(x)*ROOT::Math::chisquared_cdf(abs(x),[2]))"; }
      else                { fun="[0]*(x<[1] ? 1.0 : 0.5*ROOT::Math::chisquared_cdf_c(  x, [2]))"; }
    } else {
      if      (tsType==2) { fun=                   "[0]*ROOT::Math::chisquared_quantile_c(    x,[2])"; }
      else if (tsType==6) { fun=       "[0]*sign(0.5-x)*ROOT::Math::chisquared_quantile(sign(0.5-x)*(1.0-2.0*x),[2])"; }
      else                { fun=                   "[0]*ROOT::Math::chisquared_quantile_c(2.0*x,[2])"; }
    }
    f= new TF1 ("asym", fun, funlo, xhi);
    f->SetParameters (norm, p1, ndf, p3);
  }
  f->SetNpx(10*nbin);
  f->SetLineColor(kRed);
  return f;
}


void WorkspaceCalculator::PlotTS (TPad* canvas, TH1D* h, TH1D* ha, TH1D* hb, int tsType, Double_t obsTS, TF1** fitOut, bool norm, Color_t lineCol)
{
  if (tsType==2 && nullPOI.size()==1 && poimin[0]>=nullPOI[0]) tsType=5;   // like single-sided
  int nbin=h->GetNbinsX();
  double xlo=h->GetXaxis()->GetXmin(), xhi=h->GetXaxis()->GetXmax();
  double xw=(xhi-xlo)/nbin;
  TString tsvar=tsName(tsType), chi2label=chi2Name(tsType);

  Double_t nh= h->GetSum();

  IncludeOverflows(h);
  h->SetLineWidth(2);
  if (lineCol>=0) h->SetLineColor(lineCol);

  if (norm && nh!=0.0) h->Scale(1.0/nh);

  h->GetXaxis()->SetTitle(tsvar);
  h->GetYaxis()->SetTitle("Fraction of pseudo-experiments");
  if (tsType==5 && norm && h->GetMaximum()<0.5) h->SetMaximum (0.5 * (2.0*0.9/0.95));

  if (ha) {
    IncludeOverflows (ha);
    ha->SetLineWidth(2);
    ha->SetLineColor(kRed+1);
    Double_t nha= ha->GetSum();
    if (nha!=0.0) ha->Scale(1.0/nha);
  }

  TLine* line= 0;
  TH1* hs=0, *has=0;
  if (!isnan(obsTS) && obsTS<=xhi) {
    if (ha) {
      hs= NewShaded (h, obsTS, 3005);
      has= NewShaded (ha, obsTS, 3004);
    } else {
      hs= NewShaded (h, obsTS, -1, 7);
    }
    Double_t y= h->GetBinContent (h->FindBin(obsTS));
    if (y<0.01*h->GetMaximum()) y= 0.1*h->GetMaximum();
    line= new TLine (obsTS, 0.0, obsTS, y);
    line->SetLineWidth(3);
  }

  TF1* chi2= 0, *chi2a=0;
  if (!ha || hb) {
    chi2= AsymTF1 (tsType, 0, h, xw*(norm?1.0:nh));
    Double_t obsp= chi2->Eval(xhi);
    if (obsp>0.0) h->SetMinimum(obsp);
    if (bPOI.size()==3) {
      chi2a= dynamic_cast<TF1*>(chi2->Clone("chi2a"));
      chi2a->SetParameter (2, 2);
      chi2a->SetLineColor(kGreen);
    }
  }

  TF1* fit= 0;
  Double_t fitnorm= 1.0;
  if (hb) {
    Double_t nb= hb->GetSum();
    double xwb=(xhi-xlo)/hb->GetNbinsX();
#if 0
    // Enable this to fit to a generated chi^2
    hb->Reset();
    for (int i=0, n=int(nb); i<n; i++) {
      Double_t x= gRandom->Rndm();
      hb->Fill (x>=0.5 ? -1e-14 : ROOT::Math::chisquared_quantile(2.0*x,1.0));
    }
#endif

    Double_t fitlo= 4.0;  // Leadbetter formula only good for q0>4-6 (otherwise it is an upper bound).
    Double_t fithi= xhi;
//  Double_t fithi= isnan(obsTS) ? 0.5*xhi : 0.7*obsTS;  // use this to test performance
    // Use normalisation so fit parameters are fairly uncorrelated
    Double_t nbf= hb->Integral(hb->FindBin(fitlo),hb->FindBin(fithi));
    Double_t N= 1.0/(nb*xwb);
    Double_t Nchi= 2.0*ROOT::Math::gaussian_cdf(-sqrt(fitlo));
    fit= new TF1 ("fit", "[0]*[3]*((1.0-[1])*ROOT::Math::chisquared_pdf(x,1.0) + [4]*[1]*[2]*exp([2]*([5]-x)))", fitlo, xhi);
    fit->SetParameters (nbf, 0.0, 0.5);
    fit->SetParNames ("Nfit", "f_exp", "slope");
//  fit->FixParameter (1, 0.0);
    fit->FixParameter (2, 0.5);
    fit->FixParameter (3, 0.5/(nbf*N));
    fit->FixParameter (4, Nchi);
    fit->FixParameter (5, fitlo);
    fit->SetLineColor(kMagenta);
    if (verbose>=0) {
      cout <<"Fit "<<nbf<<" of "<<nb<<" entries in ["<<fitlo<<","<<fithi<<"], bin width "<<xwb
           <<", integral="<<fit->Integral(fitlo,fithi)/xwb<<endl;
      cout << fit->GetTitle() << " :";
      for (int i=0,n=fit->GetNpar(); i<n; i++) cout << " ["<<i<<"]="<<fit->GetParameter(i);
      cout << endl;
    }
    TFitResultPtr r = hb->Fit (fit, (verbose>=1?"LNSV":verbose>=0?"LNS":"LNSQ"), "", fitlo, fithi);
    if (verbose>=1) r->Print();
    if (verbose>=0) cout << "Integral in fit range = "<<fit->Integral(fitlo,fithi)/xwb<<endl;
    if (verbose>=0 && fit->GetNumberFreeParameters()>=2) {
      TMatrixDSym cor= r->GetCorrelationMatrix();
      if (fit->GetNumberFreeParameters()==2 && cor.GetNcols()>=2)
        cout << "Correlation = " << cor(0,1) << endl;
      else
        cor.Print();
    }
    if (!isnan(obsTS) && verbose>=0) {
      Double_t p0= N*fit->Integral(obsTS,200), ep0= N*fit->IntegralError(obsTS,200);
      double Z= -ROOT::Math::gaussian_quantile(p0,1.0);
      double eZ= max (ROOT::Math::gaussian_quantile(p0+ep0,1.0)+Z, -ROOT::Math::gaussian_quantile(p0-ep0,1.0)-Z);
      cout << Form("Fitted p-value at q0 %.2f: p0 = %.3g +/- %.3g (%.1f%%), Z = %.2f +/- %.6f",
                   obsTS, p0, ep0, ep0/p0*100.0, Z, eZ) << endl;
    }
    fitnorm= N * fit->GetParameter(3);
    fit->SetParameter (3, fitnorm * xw * (norm?1.0:nh));
  }

  Double_t xoff=0.97-canvas->GetRightMargin(), yoff=0.97-canvas->GetTopMargin();
  TLegend* leg;
  if (ha) leg= new TLegend (xoff-0.45, yoff-0.11, xoff, yoff);
  else    leg= new TLegend (xoff-0.35, yoff-0.15, xoff, yoff);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  if (0 && bPOI.size()>0 && bPOI.size()==poiName.size()) {
    TString s;
    for (size_t i=0, n=bPOI.size(); i<n; i++) {
      if (i>0) s += ", ";
      s += Form("%s=%g",poiName[i].c_str(),bPOI[i]);
    }
    leg->SetHeader (s);
  }
  if (ha) {
            leg->AddEntry (h,    Form ("%s Toys", nullModelName),       "L");
            leg->AddEntry (ha,   Form ("%s Toys",  altModelName),       "L");
  } else {
            leg->AddEntry (h,          "Pseudo-experiments",            "L");
  }
  if (line) leg->AddEntry (line,       "Observed",                      "L");
  if (chi2) leg->AddEntry (chi2, chi2label,                             "L");
  if (chi2a)leg->AddEntry (chi2a,chi2Name(tsType,2),                    "L");
  if (fit)  leg->AddEntry (fit,  Form("%s + exp fit",chi2label.Data()), "L");

#if 0
  unsigned int uoflow=0;
  if (h->GetBinContent(0)      || (ha && ha->GetBinContent(0)))      uoflow |= 1;
  if (h->GetBinContent(nbin+1) || (ha && ha->GetBinContent(nbin+1))) uoflow |= 2;
  static const char* comment[3]= {
    "Underflows in 1st bin",
    "Overflows in last bin",
    "Under/overflows in 1st/last bins"
  };
  if (uoflow) leg->AddEntry ((TObject*)0, comment[uoflow-1], "");
#endif

  if (hs && !ha && !has && hs->GetFillStyle()==1001) {
    h->Draw("axis");
    hs->Draw("same");
    h->Draw("same");
    h->Draw("axis,same");
  } else {
    h->Draw();
    if (ha) ha->Draw("same");
    if (hs) hs->Draw("same");
    if (has) has->Draw("same");
  }
  if (line) line->Draw();
  if (chi2a) chi2a->Draw ("lsame");
  if (chi2) chi2->Draw ("lsame");
  if (fit) fit->Draw("lsame");
  leg->Draw();

  TPaveText* label=0;
  if (showAtlasLabel) label= AtlasLabel (canvas, -0.70);

  PrintPlot (canvas);
  delete label;

  delete chi2;
  delete chi2a;
  delete leg;
  delete hs;
  delete line;
  if (fit) fit->SetParameter (3, fitnorm);
  if (fitOut) *fitOut= fit;
  else        delete   fit;
}


Int_t WorkspaceCalculator::ReadObject (const char* name, TObject* obj, TDirectory* dir)
{
  // Read into an already existing object.
  if (!obj || !name) return 0;
  if (!dir) dir= gROOT;
  TString msg= name;
  msg = Form ("%s object '%s'", obj->ClassName(), gSystem->PrependPathName (dir->GetPath(), msg));
  const char* dirName = gSystem->DirName(name);
  if (strcmp (dirName, ".") != 0) {
    dir= dir->GetDirectory (dirName);
    name= gSystem->BaseName(name);
  }
  Int_t ret= 0;
  if (dir)
    if (TKey* key= dir->FindKey(name))
      if (strcmp (key->GetClassName(), obj->ClassName())==0)
        ret= key->Read (obj);
  if (ret) ::Info  ("WorkspaceCalculator::ReadObject", "Read %s",           msg.Data());
  else     ::Error ("WorkspaceCalculator::ReadObject", "Could not read %s", msg.Data());
  return ret;
}


void WorkspaceCalculator::PrintResourcesUsed (const TTime& start, const char* msg)
{
  ProcInfo_t info;
  if (gSystem->GetProcInfo(&info)<0) return;
  Long_t cput= TMath::CeilNint(info.fCpuUser + info.fCpuSys);
  Long_t wall= Long64_t(gSystem->Now()-start+TTime(500))/Long64_t(1000);
  cout << Form("Resources Used:  cput=%02ld:%02ld:%02ld,mem=%ldkb,vmem=%ldkb,walltime=%02ld:%02ld:%02ld",
               cput/3600, (cput/60)%60, cput%60,
               info.fMemResident, info.fMemVirtual,
               wall/3600, (wall/60)%60, wall%60);
  if (msg) cout << " - " << msg;
  cout << endl;
}

void WorkspaceCalculator::PrintResourcesUsed(const char* msg)
{
  PrintResourcesUsed(progStart,msg);
}


bool HaveOpt (const TString& opts, char testopt, bool set=false)
{
  bool neg= false;
  for (const char* a=opts.Data(); *a; a++) {
    if      (*a=='-') neg= true;
    else if (*a=='+') neg= false;
    else if (*a==testopt) set= !neg;
  }
  return set;
}


Double_t WorkspaceCalculator::CondFit (const RooRealVar* par, Double_t val, RooAbsCollection* fitted, const char* opt)
{
  RooArgSet* poiSet= dynamic_cast<RooArgSet*>(RooArgSet(*par).snapshot());
  dynamic_cast<RooRealVar*>(poiSet->first())->setVal(val);
  Double_t nll= DoFit (*poiSet, 2, fitted, opt);
  delete poiSet;
  return nll;
}


Double_t WorkspaceCalculator::DoFit (const RooArgSet& poiSet, int type, RooAbsCollection* fitted, const char* opt)
{
  cout << "XXX do fit is run!" << endl;
  const RooArgSet* details= 0;
  Double_t nll= NaN;
  RooArgSet* poiSnap= dynamic_cast<RooArgSet*>(poiSet.snapshot());
  RooStats::ProfileLikelihoodTestStat* tss= 0;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
  RooStats::ProfileLikelihoodTestStatEnhanced* tse= dynamic_cast<RooStats::ProfileLikelihoodTestStatEnhanced*>(testStat);
  if (tse) {
    if (fitted)          tse->EnableDetailedOutput();
    ULong64_t oldopt=    tse->SetOptimize(opt,0);
    cout << "XXX spend much time enhanced!" << endl;
    nll=                 tse->EvaluateProfileLikelihood (type, *data, *poiSnap);
    if (fitted) details= tse->GetDetailedOutput();
    if (opt)             tse->SetOptimize(oldopt);
  } else
#else
  RooStats::ProfileLikelihoodTestStat* tse= 0;
#endif
  if ((tss= dynamic_cast<RooStats::ProfileLikelihoodTestStat*>(testStat))) {
    if (type==0 && HaveOpt(opt,'U',HaveOpt(optimizeTS,'U'))) type= 1;
    if (type==0 && HaveOpt(opt,'C',HaveOpt(optimizeTS,'C'))) type= 2;
    if (fitted)           tss->EnableDetailedOutput();
    cout << "XXX spend much time !" << endl;
    nll=                  tss->EvaluateProfileLikelihood (type, *data, *poiSnap);
    if (fitted)  details= tss->GetDetailedOutput();
  } else
    cerr << "No ProfileLikelihoodTestStat - skip fit" << endl;
  delete poiSnap;

  if (details) {
    TString prefix;
    if      (details->find("fitUncond_minNLL")) prefix= "fitUncond_";
    else if (details->find("fitCond_minNLL"))   prefix= "fitCond_";
    if (prefix.Length()) {
      if (tss && HaveOpt(opt,'F',HaveOpt(optimizeTS,'F'))) {
        RooArgSet fitVars;
        for (RooFIter it= details->fwdIterator(); RooAbsArg* a= it.next();) {
          const char* name= a->GetName();
          if (strncmp (name, prefix.Data(), prefix.Length()) == 0) {
            name += prefix.Length();
            a->SetName(name);  // OK to change name because we disable detailed output after
            fitVars.add(*a);
          }
        }
        ws->allVars() = fitVars;
      }
      for (RooFIter it= fitted->fwdIterator(); RooAbsArg* a= it.next();) {
        RooAbsArg* f= details->find(prefix+a->GetName());
        if (f) f->SetName(a->GetName());
      }
      *fitted = *details;
    }
  }
  if (fitted) {
    if (tse) tse->EnableDetailedOutput(0);
    if (tss) tss->EnableDetailedOutput(0);
  }
  return nll;
}


Double_t WorkspaceCalculator::UseLimits (const RooRealVar* par, Double_t val, int verbose)
{
  if        (val < par->getMin()) {
    if (verbose>=1) cout << Form("%s = %g limited by minimum at %g",par->GetName(),val,par->getMin()) << endl;
    return par->getMin();
  } else if (val > par->getMax()) {
    if (verbose>=1) cout << Form("%s = %g limited by maximum at %g",par->GetName(),val,par->getMax()) << endl;
    return par->getMax();
  }
  return val;
}


int WorkspaceCalculator::RunMinos (const RooArgSet& pars, Double_t nll_min)
{
  // If nll_min and AltMinos are specified, then assumes the
  // unconditional fit has already been done.
  if (pars.getSize() <= 0) return 0;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
  if (!(optimize & kAltMinos)) {
    if (RooStats::ProfileLikelihoodTestStatEnhanced* ts= dynamic_cast<RooStats::ProfileLikelihoodTestStatEnhanced*>(testStat)) {
      ts->UseMinos(&pars);
      Double_t nll= DoFit (RooArgSet(*pars.first()), 1, 0, "F");
      ts->UseMinos(0,kFALSE);
      return isnan(nll) ? pars.getSize() : 0;
    }
  }
#endif
  if (isnan(nll_min)) {
    nll_min= DoFit (RooArgSet(*pars.first()), 1, 0, "F");
    if (isnan(nll_min)) {
      cerr << "Unconditional fit failed - skip findSigma" << endl;
      return pars.getSize();
    }
  }

  int err= 0;
  for (RooFIter it= pars.fwdIterator(); RooAbsArg* a= it.next();) {
    RooRealVar* par= dynamic_cast<RooRealVar*>(a);
    if (!par) continue;
    err += findSigma (nll_min, par);
  }
  return err;
}



int WorkspaceCalculator::findSigma (Double_t nll_min, RooRealVar* par, Double_t nsigma,
                                    RooAbsCollection* fitted_lo, RooAbsCollection* fitted_hi,
                                    const char* opt, Double_t precision)
{
  // Find asymmetric errors for the specified parameters (pars), similar to MINOS.
  // Assumes that the parameters are at their best-fit values (for which the NLL=nll_min) and
  // the symmetric errors are a good first guess.
  RooArgSet* allVars= nullPdf->getParameters(*data);
  RooArgSet* vars= dynamic_cast<RooArgSet*>(allVars->selectByAttrib("Constant",kFALSE));
  delete allVars;
  vars->add (*par, kTRUE);  // in case it's constant
  RooArgSet* snap= dynamic_cast<RooArgSet*>(vars->snapshot());

  nsigma= fabs(nsigma);
  Double_t val= par->getVal(), err= nsigma*par->getError();
  Double_t shi= findSigma (nll_min, val+err, val, par,  nsigma, fitted_hi, opt, precision);
  *vars= *snap;
  Double_t slo= findSigma (nll_min, val-err, val, par, -nsigma, fitted_lo, opt, precision);
  *vars= *snap;
  par->setAsymError (isnan(slo)?1.0:slo, isnan(shi)?-1.0:shi);

  if (verbose>=0) {
    cout << "findSigma ";
    par->Print();
  }
  delete vars;
  delete snap;
  return (isnan(slo) ? 1 : 0) + (isnan(shi) ? 1 : 0);
}


Double_t WorkspaceCalculator::findSigma (Double_t nll_min, Double_t val_guess, Double_t val_mle, const RooRealVar* par, Double_t nsigma,
                                         RooAbsCollection* fitted, const char* opt, Double_t precision)
{
  // Find the value of sigma evaluated at a specified nsigma, assuming NLL -2logL is roughly parabolic in par.
  // The precision is specified as a fraction of the error, or based on the Minuit default tolerance.
  // Based on
  //   https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/HiggsPhys/HSG3/WWDileptonAnalysisCode/HWWStatisticsCode/trunk/macros/findSigma.C
  // by Aaron Armbruster <armbrusa@umich.edu> (2012-05-30).

  const int maxiter= 25;
  val_guess= UseLimits (par, val_guess, verbose);
  int direction = nsigma>=0.0 ? +1 : -1;
  int nDamping = 1;
  Double_t damping_factor = 1.0;
  std::map<Double_t,Double_t> guess_to_corr;
  Double_t tmu=NaN;

  if (precision <= 0.0) {
    // RooFit default tolerance is 1.0.
    // eps=0.001*tolerance. Compare with 0.005 in runPulls.C.
    Double_t eps = 0.001 * (fitTol>0.0 ? fitTol : 1.0);
    precision = 5.0*eps/(nsigma*nsigma);
  }

  int iter= 0;
  for (; iter<maxiter; iter++) {
    if (verbose>=1) cout << Form ("Parameter %s %+gsigma iteration %d: start %g (MLE%+g)", par->GetName(), nsigma, iter+1, val_guess, val_guess-val_mle) << endl;
    Double_t val_pre = val_guess;

    Double_t nll = CondFit (par, val_guess, fitted, opt);
    tmu = 2.0*(nll-nll_min);
    Double_t sigma_guess = fabs(val_guess-val_mle);
    if (tmu>0.01) sigma_guess /= sqrt(tmu);
    else          sigma_guess *= 10.0;  // protect against tmu<=0, and also don't move too far

    Double_t corr = damping_factor*(val_pre - val_mle - nsigma*sigma_guess);

    for (std::map<Double_t,Double_t>::iterator iguess= guess_to_corr.begin(); iguess != guess_to_corr.end(); iguess++) {
      if (fabs(iguess->first - val_pre) < direction*val_pre*0.02) {
        damping_factor *= 0.8;
        if (verbose>=1) cout << "Changing damping factor to " << damping_factor << endl;
        if (nDamping++ > 10) {
          nDamping = 1;
          damping_factor = 1.0;
        }
        corr *= damping_factor;
        break;
      }
    }

    // subtract off the difference in the new and damped correction
    val_guess -= corr;
    guess_to_corr[val_pre] = corr;
    val_guess= UseLimits (par, val_guess, verbose);
    Double_t relprecision = precision*fabs(val_guess-val_mle);
    Double_t delta = val_guess-val_pre;

    if (verbose>=1)
      cout << Form("%s %.3f (MLE%+.3f) -> %.3f (MLE%+.3f), change %+.3f, precision %.3f, -2lnL %.4f, sigma(guess) %.3f",
                   par->GetName(), val_pre, val_pre-val_mle, val_guess, val_guess-val_mle, delta, relprecision, tmu, sigma_guess) << endl;
    if (verbose>=2) {
      cout << "NLL:                 " << nll << endl;
      cout << "delta(NLL):          " << nll-nll_min << endl;
      cout << "nsigma*sigma(pre):   " << fabs(val_pre-val_mle) << endl;
      cout << "sigma(guess):        " << sigma_guess << endl;
      cout << "par(guess):          " << val_guess+corr << endl;
      cout << "best-fit val:        " << val_mle << endl;
      cout << "tmu:                 " << tmu << endl;
      cout << "Precision:           " << direction*val_guess*precision << endl;
      cout << "Correction:          " << -corr << endl;
      cout << "nsigma*sigma(guess): " << fabs(val_guess-val_mle) << endl;
      cout << endl;
    }

    if (fabs(delta) <= relprecision) break;
  }
  if (iter >= maxiter) {
    cerr << "findSigma failed after " << iter << " iterations" << endl;
    return NaN;
  }

  Double_t err= val_guess-val_mle;
  if (verbose>=0) cout << Form ("%s %+gsigma = %.3f at -2lnL = %.4f after %d iterations", par->GetName(), nsigma, err, tmu, iter+1) << endl;

  return err;
}


//==============================================================================
// Constructors and destructor
//==============================================================================

WorkspaceCalculator::WorkspaceCalculator (const char* name)
  : TNamed (name, name)
{
  Init();
}

WorkspaceCalculator::~WorkspaceCalculator()
{
}

int WorkspaceCalculator::ParseArgs (const char* args)
{
  char* str= strcpy (new char[strlen(args)+1], args);
  const char* argv[4096]= {GetName()};
  int argc= 0, maxarg= sizeof(argv)/sizeof(argv[0])-1;
  for (char* a=str; (argv[++argc]= strtok (a, " ")) && argc<maxarg; a= NULL) ;
  int ret= ParseArgs (argc, argv);
  delete [] str;
  return ret;
}
