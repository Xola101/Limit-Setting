/*==============================================================================
$Id: LikelihoodCalculator.cxx 697232 2015-09-28 19:37:09Z adye $

Top-level class to perform a parameter scan to obtain the likelihood curve or contours.
This calculator is used by LikelihoodScan.

Main parameters:

  fileName         workspace file or list of result files
  wsName           workspace name
  modelSBName      ModelConfig object name
  dataName         observed dataset

  plotResult:      plot result of test
  writeResult:     write result (default is true)

Author: Tim Adye.

==============================================================================*/

#include "Utilities/WorkspaceCalculatorConfig.h"
#include "LikelihoodCalculator.h"

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
#include "TList.h"
#include "TObjArray.h"
#include "TKey.h"
#include "RooRandom.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
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
#include "TMarker.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStopwatch.h"

#ifdef USE_THn
#include "THn.h"
#else
#include "THnSparse.h"
#endif

#include "Math/DistFunc.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"


using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::min;
using std::max;
using std::isnan;

std::vector<Double_t> LikelihoodCalculator::GetUpLevels (int nupLevels, int ndim, int verbose, bool lim95)
{
  vector<Double_t> levels(nupLevels);
  for (int sig=1; sig<=nupLevels; sig++) {
    Double_t up;
    if (lim95 && sig==2) {
      up= ROOT::Math::chisquared_quantile (0.95, ndim);
      if (verbose>=1) cout << Form("%dD p-value %.11g = UP %g", ndim, 0.95, up) << endl;
    } else if (ndim==1) {   // handle simple case more accurately
      up= sig*sig;
    } else {
      double psig= 2.0*ROOT::Math::gaussian_cdf(sig)-1.0;
      up= ROOT::Math::chisquared_quantile (psig, ndim);
      if (verbose>=1) cout << Form("%dD %d sigma = p-value %.11g = UP %g", ndim, sig, psig, up) << endl;
    }
    levels[sig-1]= up;
  }
  return levels;
}

const char* percent (Double_t p, bool sign=false)
{
  return Form (Form ("%%%s.%df%%%%", sign?"+":"", max(0,-int(floor(log10(1.0-fabs(p))+2.5)))), 100.0*p);
}

const char* UpLevelName (Double_t up, int ndim=1, int style=0, Double_t scale=0.65, const char* latex="#", bool sign=false)
{
  Double_t Z, p;
  if (ndim==1) {
    Z= sqrt(fabs(up));
    p= 1.0-2.0*ROOT::Math::gaussian_cdf(-Z);
  } else {
    p= ROOT::Math::chisquared_cdf (fabs(up), double(ndim));
    Z= -ROOT::Math::gaussian_quantile (0.5*(1.0-p), 1.0);
  }
  Int_t sig= TMath::Nint(Z);
  if (style ? (style==1) : (Z>4.0 || fabs(Z-sig) < 1e-6)) {
    if (sig==0) return ("MLE");
    if (up<0.0) sig= -sig;
    return Form(Form("%%%sd%%ssigma",sign?"+":""), sig, latex);
  } else {
    if (up<0.0) p= -p;
    const char* s= percent (p, sign);
    if (!scale || strcmp(latex,"#")!=0) return s;
    return Form("#scale[%g]{%s}",scale,s);
  }
}

void LikelihoodCalculator::Init() {
  optimize           &= ~kSkipDataTS;
  seed               = -1;    // don't use random numbers, so no point in randomising and saving seeds
  resultFileName     = "scan.root";
  resultName         = "scan";

  jobSet             = 0;
  nscanDefault       = 40;    // Default number of mu points to scan in [scanMin:scanMax]
  nscan              = nscanDefault;
  npoi               = 1;
  firstPoint         = 0;
  lastPoint          = -1;    // Set >=0 to scan only a subset of the points
  runScan            = true;

  res                = 0;
  interpolate        = 1;
  plotAllProjections = 0;
  ymax               = 7.0;
  ncontour2D         = 2;
  nupLevels          = 3;
  lim95              = false;
  massAxis           = -1;
  plotOpt            = "pl";
  toPlot.SetOwner();
}


void LikelihoodCalculator::SetDefaults() {
  if (testStatType<0) testStatType= 2; // test statistic type (2=two-sided)
  WorkspaceCalculator::SetDefaults();
}


void LikelihoodCalculator::ShowParms() const
{
     cout << GetName()
          << " wsName=\""         << wsName      << "\""
          << ", modelSBName=\""   << modelSBName << "\""
          << ", modelBName=\""    << modelBName  << "\""
          << ", dataName=\""      << dataName    << "\"";
     PrintOptimize (optimize & ~kAdjustRanges, ", ");
     if (optimize & kAdjustRanges)   cout << ", AdjustRanges "<<newParameterRanges<<" sigma";
     if (detailedOutput)             cout << ", DetailedOutput";
     if (dropNegative)               cout << ", dropNegative";
     cout << ", npoints="         << nscan;
     if (npoi!=1) cout << " (" << Join (npoints, "x") << ")";
     int nv=0;
     for (size_t d=0, nd=pointsToScan.size(); d<nd; d++) nv += pointsToScan[d].size();
     if (nv>0) {
       cout << ", pointsToScan=";
       for (size_t d=0, nd=pointsToScan.size(); d<nd; d++)
         cout << (d>0 ? " : " : "") << Join(pointsToScan[d]);
     }
     if      (lastPoint==firstPoint)
     cout << ", scanPoint="       << firstPoint;
     else if (lastPoint> firstPoint)
     cout << ", scanPoints="      << firstPoint << '-' << lastPoint;
     if (invMass>0.0)
     cout << ", invMass="         << invMass;
     if (!poimin.empty())
     cout << ", poimin="          << Join(poimin);
     if (!poimax.empty())
     cout << ", poimax="          << Join(poimax);
     if (!scanMin.empty())
     cout << ", scanMin="         << Join(scanMin);
     if (!scanMax.empty())
     cout << ", scanMax="         << Join(scanMax);
     cout << endl;
}


int LikelihoodCalculator::RunOverWorkspace()
{
  // run the calculator
  int     ok= SetupWorkspace();
  if (ok) ok= SetupMinimizer();
  if (ok) ok= SetupModel();
  if (ok) ok= SetupInitialParms();
  if (ok) ok= SetupTestStat();
  if (ok) ok= RunCalculatorProof();
  if (ok) ok= GetTestStatInfo();
  DeleteCalculator();   // don't use from here on
  return ok;
}


TAxis* LikelihoodCalculator::ReplaceAxis (TAxis* dest, const TAxis* source)
{
  // Replaces an axis with that of a different histogram
  TString  name= dest->GetName();
  TObject* hist= dest->GetParent();
  source->Copy (*dest);
  dest->SetName   (name);
  dest->SetParent (hist);
  return dest;
}

vector<Double_t> LikelihoodCalculator::GetBinPoints (const TAxis* ax)
{
  // calculate scan points from bin edges, assuming edges chosen at mid-point between points.
  Int_t nb= ax->GetNbins();
  vector<Double_t> points(nb);
  if (ax->IsVariableBinSize()) {
    const Double_t* b= ax->GetXbins()->GetArray();
    points[0]= 0.5*(b[0]+b[1]);
    for (Int_t i=1; i<nb; i++) {
      points[i]= 2.0*b[i] - points[i-1];
      if (points[i] <= b[i] || points[i] >= b[i+1]) points[i]= 0.5*(b[i]+b[i+1]);  // stay sane
    }
  } else {
    for (Int_t i=0; i<nb; i++) points[i]= ax->GetBinCenter(i+1);
  }
  return points;
}

Double_t LikelihoodCalculator::GetBinPoint (const TAxis* ax, Int_t bin)
{
  // calculate scan points from bin edges, assuming edges chosen at mid-point between points.
  if (ax->IsVariableBinSize() && bin >= 1 && bin <= ax->GetNbins()) {
    const Double_t* b= ax->GetXbins()->GetArray();
    Double_t pt= 0.5*(b[0]+b[1]);
    for (Int_t i=1; i<bin; i++) {
      pt= 2.0*b[i] - pt;
      if (pt <= b[i] || pt >= b[i+1]) pt= 0.5*(b[i]+b[i+1]);  // stay sane
    }
    return pt;
  } else
    return ax->GetBinCenter(bin);
}

Double_t LikelihoodCalculator::GetBinPoint (const HistTypeBase* h, Int_t d, Int_t bin)
{
  return GetBinPoint (h->GetAxis(d), bin);
}

std::vector<int>
LikelihoodCalculator::GetBin (int ibin, const std::vector<int>& npoints, bool lr, int offset, int dfix, int vfix)
{
  int npoi=npoints.size();
  vector<int> ipt(npoi);
  for (int d=0; d<npoi; d++) {
    if (d==dfix) { ipt[d]= vfix; continue; }
    int np= npoints[d], id= ibin % np;
    ibin /= np;
    if (lr && ibin%2) id= np-id-1;  // scan left-right-left-right-... for faster fits
    ipt[d]= id+offset;
  }
  return ipt;
}


int
LikelihoodCalculator::GetBin (const std::vector<int>& bin, const std::vector<int>& npoints, bool lr, int offset, int dfix)
{
  int npoi=npoints.size(), ibin=0, span=1;
  for (int d=0; d<npoi; d++)
    if (d!=dfix) span *= npoints[d];
  for (int d=npoi-1; d>=0; d--) {
    if (d==dfix) continue;
    int np= npoints[d], id= bin[d] - offset;
    if      (id<0)   id= 0;
    else if (id>=np) id=np-1;
    if (lr && ibin%2) id= np-id-1;  // scan left-right-left-right-... for faster fits
    span /= np;
    ibin += span*id;
  }
  return ibin;
}


TString LikelihoodCalculator::BinName (const HistTypeBase* h, const std::vector<Int_t>& bin, Int_t skip)
{
  // Name and location for histogram global bin number
  Int_t npoi= h->GetNdimensions();
  TString st, sp;
  for (int d=0; d<npoi; d++) {
    if (d==skip) continue;
    if (d>0) st += ",";
    st += h->GetAxis(d)->GetTitle();
    if (d>0) sp += ",";
    sp += Form("%g",GetBinPoint(h,d,bin[d]));
  }
  return st + " = " + sp;
}


TString LikelihoodCalculator::HistShape (const HistTypeBase* h)
{
  // Show histogram shape
  Int_t npoi= h->GetNdimensions();
  TString st;
  Long64_t nb= 1;
  for (int d=0; d<npoi; d++) {
    if (d>0) st += ",";
    const TAxis* ax= h->GetAxis(d);
    st += Form ("%dx%g:%g",ax->GetNbins(),ax->GetXmin(),ax->GetXmax());
    nb *= ax->GetNbins();
  }
  st += Form("(%g/%lld)",h->GetEntries(),nb);
  return st;
}


int LikelihoodCalculator::RunCalculator()
{
   // scan likelihood for different values of the POIs.
   TStopwatch tw;

   RooArgSet* poiSet= dynamic_cast<RooArgSet*>(nullModel->GetParametersOfInterest()->snapshot());
   RooArgList pois= *poiSet;

   if (npoi>pois.getSize()) npoi= pois.getSize();
   if (npoints.empty()) npoints.resize(1,nscanDefault);
   npoints.resize(npoi,1);
   fitInRange.resize(npoi);
   res2.resize(npoi);
   pointsToScan.resize(npoi);

   vector<Double_t> step(npoi);
   scanMin.resize(npoi,NaN);
   scanMax.resize(npoi,NaN);
   nscan= 1;   // recalculate
   for (size_t i= 0; i<size_t(npoi); i++) {
     nscan *= npoints[i];
     RooRealVar* poi = dynamic_cast<RooRealVar*>(pois.at(i));
     if (isnan(scanMin[i])) scanMin[i]= poi->getMin();
     if (isnan(scanMax[i])) scanMax[i]= poi->getMax();
     if (!pointsToScan[i].empty()) continue;
     step[i]= (scanMax[i]-scanMin[i])/npoints[i];
     if (fitInRange[i]) {
       poiSet->remove(*poi);  // don't fix for conditional fit
       RooArgSet* vars= nullPdf->getParameters(*data);
       RooAbsArg* var= vars->find(*poi);
       delete vars;
       if (!var) {
         cout << "Variable "<<poi->GetName()<<" not found in PDF"<<endl;
         delete poiSet;
         return 0;
       }
       pois.replace (*poi, *var);
     }
     if (verbose>=0)
       cout << "Scan "<<poi->GetName()<<" from "<<scanMin[i]<<" to "<<scanMax[i]<<" in steps of "<<step[i]<<endl;
   }

   if (firstPoint>=nscan) {
     delete poiSet;
     return 0;
   }
   if (lastPoint>=nscan || lastPoint<firstPoint) lastPoint = nscan-1;

   RooArgSet xset;
   if (size_t n=extraPOI.size()) {
     RooArgSet allVars= ws->allVars();
     for (size_t i=0; i<n; i++) {
       RooAbsCollection* vars= allVars.selectByName(extraPOI[i].c_str());
       if (vars->getSize())
         xset.add(*vars);
       else
         cerr << "Extra POI '"<<extraPOI[i]<<"' not found" << endl;
       delete vars;
     }
     if (verbose>=0) {
       cout << "The following variables will follow the POI scan:-" << endl;
       xset.Print("v");
     }
   }

   res= new HistType (resultName, "-lnL", npoi, &npoints[0], &scanMin[0], &scanMax[0]);
   for (int d=0; d<npoi; d++) {
     res->GetAxis(d)->SetTitle(pois.at(d)->GetName());
     if (!pointsToScan[d].empty()) {
       size_t n= pointsToScan[d].size();
       vector<Double_t> lo(n+1);
       for (size_t i=1; i<n; i++)
         lo[i]= 0.5*(pointsToScan[d][i-1] + pointsToScan[d][i]);
       lo[0]= 2.0*pointsToScan[d][0]   - lo[1];
       lo[n]= 2.0*pointsToScan[d][n-1] - lo[n-1];
       res->GetAxis(d)->Set (n, &lo[0]);
       res->GetAxis(d)->SetLimits(scanMin[d],scanMax[d]);  // put back limits, which may not be the same as 1st/last bin
     }
     if (!fitInRange[d]) continue;
     res2[d]= new HistType (Form("%s_%s",resultName.Data(),pois.at(d)->GetName()), pois.at(d)->GetName(),
                            npoi, &npoints[0], &scanMin[0], &scanMax[0]);
   }
   for (int d=0; d<npoi; d++) {
     if (res2[d])
       for (int d2=0; d2<npoi; d2++)
         ReplaceAxis (res2[d]->GetAxis(d2), res->GetAxis(d2));
   }

   for (int i= firstPoint; i<=lastPoint; i++) {
     vector<Double_t> pt(npoi);
     vector<int> ipt= GetBin (i, npoints, true, 0);
     for (int d= 0; d<npoi; d++) {
       if (!pointsToScan[d].empty()) pt[d]= pointsToScan[d][ipt[d]];
       else                          pt[d]= scanMin[d] + (ipt[d]+0.5)*step[d];
       RooRealVar* poi= dynamic_cast<RooRealVar*>(pois.at(d));
       if (!fitInRange[d]) poi->setVal(pt[d]);
       else {
         Double_t lo, hi;
         if (!pointsToScan[d].empty()) {
           lo= ipt[d]>0                             ? 0.5 * (pt[d] + pointsToScan[d][ipt[d]-1]) : poimin[d];
           hi= ipt[d]<int(pointsToScan[d].size()-1) ? 0.5 * (pt[d] + pointsToScan[d][ipt[d]+1]) : poimax[d];
         } else {
           lo= pt[d] - 0.5*step[d];
           hi= pt[d] + 0.5*step[d];
         }
         if (i==firstPoint || lo != poi->getMin() || hi != poi->getMax()) {  // otherwise leave value from last fit
           poi->setRange (lo, hi);
           if (poi->getError() > 0.45*(hi-lo)) poi->setError(0.45*(hi-lo));
           if      (d<int( bPOI.size()) &&  bPOI[d] >= lo &&  bPOI[d] <= hi) poi->setVal( bPOI[d]);
           else if (d<int(sbPOI.size()) && sbPOI[d] >= lo && sbPOI[d] <= hi) poi->setVal(sbPOI[d]);
           else poi->setVal(pt[d]);
         }
       }
     }
     if (verbose>=1 || (verbose>=0 && i==firstPoint)) {
       if (verbose==0) cout << "First of "<<lastPoint-firstPoint+1<<" scan points: ";
       else            cout << "Scan point ";
       cout <<i<<" at "<<pois.at(0)->GetName();
       for (int d= 1; d<npoi; d++) cout << ',' << pois.at(d)->GetName();
       cout <<" = "<<ipt[0];
       for (int d= 1; d<npoi; d++) cout << ',' << ipt[d];
       cout << " = ";
       for (int d= 0; d<npoi; d++) {
         if (d>0) cout << ',';
         if (fitInRange[d]) {
           RooRealVar* poi= dynamic_cast<RooRealVar*>(pois.at(d));
           cout << poi->getMin() << ':' << poi->getMax();
         } else
           cout << pt[d];
       }
       cout << endl;
       pois.Print("v");
     }
     for (RooLinkedListIter it= xset.iterator(); TObject* o= it.Next();)
       if (RooRealVar* v= dynamic_cast<RooRealVar*>(o))
         v->setVal(pt[0]);   // extra POIs only follow first POI
     if (verbose>=1) {
       cout << "Set variables to the POI value:-" << endl;
       xset.Print("v");
     }

     // use type=2 to skip unconditional fit
     Double_t nll= DoFit (*poiSet, 2);
     vector<int> ibin(npoi);
     for (int d=0;d<npoi;d++) ibin[d]= ipt[d]+1;
     res->SetBinContent(&ibin[0],nll);
     for (int d=0; d<npoi; d++) {
       if (res2[d]) res2[d]->SetBinContent (&ibin[0], dynamic_cast<RooRealVar*>(pois.at(d))->getVal());
     }
   }

   // restore limits for fitInRange POIs
   for (int d= 0; d<npoi; d++) {
     if (!fitInRange[d]) continue;
     RooRealVar* poi= dynamic_cast<RooRealVar*>(pois.at(d));
     poi->setMin (poimin[d]);
     poi->setMax (poimax[d]);
   }

   delete poiSet;

   if (verbose>=0) {
     cout << "Time to scan likelihood function: "; tw.Print();
   }
   return 1;
}


Int_t LikelihoodCalculator::FindBinExt (const TAxis* ax, Double_t x, bool upper)
{
  // Like TAxis::FindFixBin, but allows bins outside 1..N and handle variable bins where min/max don't correspond to 1st/last bin
  Int_t nb= ax->GetNbins();
  Double_t xlo= ax->GetBinLowEdge(1), xhi= ax->GetBinUpEdge(nb);
  Int_t bin;
  if        (x<xlo)
    bin= TMath::FloorNint((x-xlo)/ax->GetBinWidth(1))  + 1;
  else if (!(x<xhi)) {
    bin= TMath::FloorNint((x-xhi)/ax->GetBinWidth(nb)) + 1 + nb;
  } else
    bin= ax->FindFixBin(x);
  if (upper && x==GetBinLowEdgeExt(ax,bin)) bin--;
  return bin;
}


Double_t LikelihoodCalculator::GetBinLowEdgeExt (const TAxis* ax, Int_t bin)
{
  // Like TAxis::GetBinLowEdge, but allows bins outside 1..N and handle variable bins where min/max don't correspond to 1st/last bin
  if (bin<1)  return ax->GetBinLowEdge(1) + (bin-1)   *ax->GetBinWidth(1);
  Int_t nb= ax->GetNbins();
  if (bin>nb) return ax->GetBinUpEdge(nb) + (bin-1-nb)*ax->GetBinWidth(nb);
  return ax->GetBinLowEdge(bin);
}


int LikelihoodCalculator::AddMassAxis (HistTypeBase*& r, std::vector<HistTypeBase*>& r2) const
{
  Int_t nd0=r->GetNdimensions(), nd=nd0+1;
  size_t nv= min (min(poiName.size(), npoints.size()), min(scanMin.size(), scanMax.size()));
  Int_t dm=-1;
  for (size_t d=0; d<nv; d++)
    if (poiName[d]=="M" && npoints[d]>0 && !isnan(scanMin[d]) && !isnan(scanMax[d])) {
      dm=d; break;
    }
  if (dm<0) return dm;

  TAxis* axm;
  if (pointsToScan[dm].empty())
    axm= new TAxis(npoints[dm],scanMin[dm],scanMax[dm]);
  else {
    size_t n= pointsToScan[dm].size();
    vector<Double_t> lo(n+1);
    for (size_t i=1; i<n; i++)
      lo[i]= 0.5*(pointsToScan[dm][i-1] + pointsToScan[dm][i]);
    lo[0]= 2.0*pointsToScan[dm][0]   - lo[1];
    lo[n]= 2.0*pointsToScan[dm][n-1] - lo[n-1];
    axm= new TAxis (n, &lo[0]);
    axm->SetLimits(scanMin[dm],scanMax[dm]);  // put back limits, which may not be the same as 1st/last bin
  }
  axm->SetName(Form("axis%d",dm));
  axm->SetTitle(poiName[dm].c_str());
  Int_t massBin= axm->FindFixBin(invMass);
  if (verbose>=0)
    cout << Form("New axis %d with bins %dx%g:%g for mass. Fill bin %d with mass %g GeV",
                 dm, npoints[dm], scanMin[dm], scanMax[dm], massBin, invMass) << endl;

  vector<Double_t> lo(nd), hi(nd);
  vector<Int_t> nb(nd), id(nd);
  for (Int_t d=0,e=0; d<nd; d++) {
    const TAxis* ax= (d==dm) ? axm : r->GetAxis(e);
    id[d]= (d==dm ? -1 : e++);
    nb[d]= ax->GetNbins();
    lo[d]= ax->GetXmin();
    hi[d]= ax->GetXmax();
  }

  HistType* h= new HistType (r->GetName(), r->GetTitle(), nd, &nb[0], &lo[0], &hi[0]);
  vector<HistTypeBase*> h2(nd);
  for (Int_t d=0; d<nd; d++) {
    Int_t e=id[d];
    ReplaceAxis (h->GetAxis(d), (e<0 ? axm : r->GetAxis(e)));
    if (e<0 || !r2[e]) continue;
    h2[d]= new HistType (r2[e]->GetName(), r2[e]->GetTitle(), nd, &nb[0], &lo[0], &hi[0]);
    for (int d2=0; d2<nd; d2++) {
      Int_t e2=id[d2];
      ReplaceAxis (h2[d]->GetAxis(d2), (e2<0 ? axm : r2[e]->GetAxis(e2)));
    }
  }

  for (int i= 0, n=r->GetNbins(); i<n; i++) {
    vector<Int_t> bold(nd0), bnew(nd);
    Double_t v= r->GetBinContent(i,&bold[0]);
    for (Int_t d=0; d<nd; d++) {
      Int_t e=id[d];
      bnew[d]= (e<0) ? massBin : bold[e];
    }
    h->SetBinContent(&bnew[0],v);
    for (Int_t d=0; d<nd; d++) {
      int e=id[d];
      if (e>=0 && r2[e]) h2[d]->SetBinContent (&bnew[0], r2[e]->GetBinContent(&bold[0]));
    }
  }

  delete r;
  for (Int_t d=0; d<nd0; d++) delete r2[d];
  r= h;
  r2= h2;
  return dm;
}


int LikelihoodCalculator::ResizeResult (HistTypeBase*& r, std::vector<HistTypeBase*>& r2) const
{
  // Recreate histogram with a new range, but keep same bin widths

  if (scanMin.empty() && scanMax.empty()) return 0;

  Int_t nd=r->GetNdimensions();
  vector<Double_t> lo=scanMin, hi=scanMax;
  lo.resize(nd,NaN);
  hi.resize(nd,NaN);
  vector<Int_t> nb(nd), off(nd), mod(nd);
  int nmod= 0;
  for (Int_t d=0; d<nd; d++) {
    const TAxis* ax= r->GetAxis(d);
    if (isnan(lo[d])) lo[d]= ax->GetXmin();
    else if (ax->GetXmin() != lo[d]) mod[d]++;
    if (isnan(hi[d])) hi[d]= ax->GetXmax();
    else if (ax->GetXmax() != hi[d]) mod[d]++;
    if (!mod[d]) continue;
    Int_t ilo= FindBinExt(ax,lo[d]);
    Int_t ihi= FindBinExt(ax,hi[d],true);
    off[d]= ilo-1;
    nb[d]= ihi-ilo+1;
    // adjust lo/hi to bin edges
    lo[d]= GetBinLowEdgeExt(ax,ilo);
    hi[d]= GetBinLowEdgeExt(ax,ihi+1);
    if (verbose>=0)
      cout << Form("Adjust %s axis %s bins %dx%g:%g -> %dx%g:%g (%d:%d)",
                   ax->GetTitle(), (ax->IsVariableBinSize() ? "variable" : "fixed"),
                   ax->GetNbins(), ax->GetXmin(), ax->GetXmax(), nb[d], lo[d], hi[d], ilo, ihi) << endl;
    nmod++;
  }
  if (!nmod) return 0;

  HistType* h= new HistType (r->GetName(), r->GetTitle(), nd, &nb[0], &lo[0], &hi[0]);
  for (Int_t d=0; d<nd; d++) {
    const TAxis* ar= r->GetAxis(d);
          TAxis* ah= h->GetAxis(d);
    if (!mod[d])
      ReplaceAxis (ah, ar);
    else {
      if (ar->IsVariableBinSize()) {
        if (off[d]<0 || off[d]+nb[d]>ar->GetNbins()) {
          Int_t n= nb[d]+1;
          vector<Double_t> bin(n);
          for (Int_t i=0; i<n; i++)
            bin[i]= GetBinLowEdgeExt(ar,i+off[d]+1);
          ah->Set (nb[d], &bin[0]);
        } else
          ah->Set (nb[d], ar->GetXbins()->GetArray() + off[d]);
        ah->SetLimits(lo[d],hi[d]);
      }
      ah->SetTitle (r->GetAxis(d)->GetTitle());
    }
  }

  vector<HistTypeBase*> h2(nd);
  for (int d=0; d<nd; d++) {
    if (!r2[d]) continue;
    h2[d]= new HistType (r2[d]->GetName(), r2[d]->GetTitle(), nd, &nb[0], &lo[0], &hi[0]);
    for (int d2=0; d2<nd; d2++)
      ReplaceAxis (h2[d]->GetAxis(d2), h->GetAxis(d2));
  }

  for (Int_t i= 0, n=r->GetNbins(); i<n; i++) {
    vector<Int_t> bold(nd), bnew(nd);
    Double_t v= r->GetBinContent(i,&bold[0]);
    Int_t d=0;
    for (; d<nd; d++) {
      bnew[d]= bold[d]-off[d];
      if (bnew[d]<1 || bnew[d]>nb[d]) break;
    }
    if (d<nd) continue;
    h->SetBinContent(&bnew[0],v);
    for (d=0; d<nd; d++) {
      if (r2[d]) h2[d]->SetBinContent (&bnew[0], r2[d]->GetBinContent(&bold[0]));
    }
  }
  delete r;
  for (Int_t d=0; d<nd; d++) delete r2[d];
  r= h;
  r2= h2;
  return nmod;
}


void LikelihoodCalculator::AdjustResults (bool resultWasRead)
{
  if (!res) return;
  if (!resultWasRead) {
    saved.Add (res);
    for (int d=0, n=res2.size(); d<n; d++)
      if (res2[d]) saved.Add(res2[d]);
    return;
  }

  npoi= res->GetNdimensions();
  poiName.resize(npoi);
  npoints.resize(npoi);
  fitInRange.resize(npoi);
  scanMin.resize(npoi);
  scanMax.resize(npoi);
  pointsToScan.resize(npoi);

  if (massAxis>=0) {
    poiName[massAxis]= "m_{H} [GeV]";
    invMass= 0.0;
  }
  nscan= 1;
  for (int d=0; d<npoi; d++) {
    TAxis* ax= res->GetAxis(d);
    const char* title= ax->GetTitle();
    if       (!poiName[d].size()) {
      poiName[d]= title;
    } else if (poiName[d] != title) {
      string oldName= title;
      cout << "Rename " << title << " -> " << poiName[d] << endl;
      res->GetAxis(d)->SetTitle(poiName[d].c_str());
      if (res2[d]) res2[d]->SetTitle (poiName[d].c_str());
      for (int d2=0; d2<npoi; d2++)
        if (res2[d2]) res2[d2]->GetAxis(d)->SetTitle(poiName[d].c_str());
      poiName[d]= oldName;  // poiName holds original name from now on
    }

    npoints[d]= ax->GetNbins();
    scanMin[d]= ax->GetXmin();
    scanMax[d]= ax->GetXmax();
    nscan *= npoints[d];
    if (ax->IsVariableBinSize()) {
      pointsToScan[d]= GetBinPoints (ax);
      const Double_t* b= ax->GetXbins()->GetArray();
      cout << Form("%s %dx%g:%g bins: ", ax->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax())
           << Join(pointsToScan[d],",",false)
           << "\n  bin edges: " << Join(vector<Double_t>(b,b+ax->GetXbins()->GetSize()),",",false) << endl;
    }
  }

  saved.Add (res);
  for (int d=0; d<npoi; d++)
    if (res2[d]) saved.Add(res2[d]);
}

//==============================================================================
// Other methods
//==============================================================================

TTree* LikelihoodCalculator::GetLimits()
{
  int nPoints=0, nFail=0;

  UInt_t initialSeed= seed, fitInRangeBits=0;
  for (size_t i=0, n=fitInRange.size(); i<n; i++)
    if (fitInRange[i]) fitInRangeBits |= (1<<i);
  TTree* bandTree= new TTree ("band", GetName());
  bandTree->SetDirectory(0);   // memory resident ntuple
  bandTree->Branch ("invMass",      &invMass);
  bandTree->Branch ("nToys",        &nPoints);
  bandTree->Branch ("nFail",        &nFail);
  bandTree->Branch ("muhat",        &poihat);
  bandTree->Branch ("muhat_err",    &poierr);
  bandTree->Branch ("muMin",        &scanMin);
  bandTree->Branch ("muMax",        &scanMax);
  bandTree->Branch ("tstype",       &testStatType);
  bandTree->Branch ("wsfile",       &wsfile);
  bandTree->Branch ("optimize",     &optimize);
  bandTree->Branch ("seed",         &initialSeed);
  bandTree->Branch ("fitInRange",   &fitInRangeBits);
  bandTree->Branch ("args",         &cmdArgs);
  bandTree->Fill();   // just one entry so we can combine later
  bandTree->ResetBranchAddresses();

  saved.Add (bandTree);
  return bandTree;
}


int LikelihoodCalculator::ReadResult (TFile* f, const char* objname)
{
  HistTypeBase* r= ReadResult (f, objname, verbose);
  if (!r) return 1;
  int ndim= r->GetNdimensions();
  vector<HistTypeBase*> r2(ndim);
  for (int d=0; d<ndim; d++) {
    TString objname2= Form("%s_%s",objname,r->GetAxis(d)->GetTitle());
    if (!f->FindKey(objname2)) continue;
    if (verbose>=0) cout << "Read '"<<objname2<<"' from file " << f->GetName() << endl;
    r2[d]= ReadResult (f, objname2, verbose);
  }

  massAxis= AddMassAxis(r,r2);
  int nmod= ResizeResult(r,r2);
  ndim= r->GetNdimensions();  // could have been bumped by AddMassAxis

  int err= MergeResult (res, r, f->GetName());
  if (verbose>=1) {
    if (massAxis || nmod) cout << "Adjusted " << HistShape(r)   << endl;
                          cout << "Merged   " << HistShape(res) << endl;
  }
  delete r;
  if (int(res2.size())<ndim) res2.resize(ndim);
  for (int d=0; d<ndim; d++) {
    if (!r2[d]) continue;
    MergeResult (res2[d], r2[d], f->GetName());
    delete r2[d];
  }
  return err;
}


LikelihoodCalculator::HistTypeBase*
LikelihoodCalculator::ReadResult (TFile* f, const char* objname, int verbose)
{
  HistTypeBase* radd = 0;
  f->GetObject (objname, radd);
  if (!radd) {
    TH1* hadd = 0;  // Check for a TH1, which is what we used to use.
    f->GetObject (objname, hadd);
    if (!hadd) return 0;
#ifdef USE_THn
    radd= HistType::CreateHn     (objname, hadd->GetTitle(), hadd);
#else
    radd= HistType::CreateSparse (objname, hadd->GetTitle(), hadd);
#endif
    delete hadd;
    if (!radd) return 0;
  }
  if (verbose>=1)
    cout << radd->ClassName() << "::" << radd->GetName() << " is " << HistShape(radd) << endl;
  return radd;
}


int LikelihoodCalculator::MergeResult (HistTypeBase*& res, const HistTypeBase* radd, const char* fname)
{
  int err=0;
  if (!res) {
    res= dynamic_cast<HistType*>(radd->Clone());
  } else {
    Int_t ndim= res->GetNdimensions();
    for (int d=0; d<ndim; d++) {
      if (radd->GetNdimensions()       != ndim                        ||
          radd->GetAxis(d)->GetNbins() != res->GetAxis(d)->GetNbins() ||
          radd->GetAxis(d)->GetXmin()  != res->GetAxis(d)->GetXmin()  ||
          radd->GetAxis(d)->GetXmax()  != res->GetAxis(d)->GetXmax()) {
        cerr << "Histogram shape is different in file " << fname
             << ": " << HistShape(res) << " -> " << HistShape(radd) << endl;
        return 2;
      }
      if (strcmp (radd->GetAxis(d)->GetTitle(), res->GetAxis(d)->GetTitle()) != 0) {
        cerr << "POI name changed from " << res->GetAxis(d)->GetTitle()
             << " to " << radd->GetAxis(d)->GetTitle() << " in file " << fname<<endl;
        err= 2;
      }
    }
    for (int i= 0, n=radd->GetNbins(); i<n; i++) {
      vector<Int_t> bin(ndim);
      Double_t vnew= radd->GetBinContent(i,&bin[0]);
      if (vnew==0.0) continue;
      Double_t vold= res ->GetBinContent(&bin[0]);
      if (vold!=0.0 && fabs(vnew-vold)>1e-4) {
        TString binName= BinName(res,bin);
        cerr << Form("Point %d at %s changed from %.4f to %.4f in file %s",
                     i, binName.Data(), vold, vnew, fname)<<endl;
        err=2;
      }
      res->SetBinContent(&bin[0],vnew);
    }
  }
  return err;
}


int LikelihoodCalculator::ReadTree (TFile* f, UInt_t& thisSeed)
{
  thisSeed= 0;
  TTree* bandTree = 0;
  f->GetObject ("band", bandTree);
  if (!bandTree) return 0;
  Int_t nent = bandTree->GetEntries();
  if (nent != 1) {
    cerr << "results file contains " << nent << " entries" << endl;
    if (!force) return 1;
    if (nent==0) return 0;
  }

  Double_t invMass1=0.0, poihat1=0.0, poierr1=0.0;
  UInt_t fitInRangeBits=(UInt_t)-1;
  Int_t testStatType1=-1;
  string* wsfile1=0;
  ULong64_t optimize1=0;
  bandTree->SetBranchAddress ("invMass", &invMass1);
  if (bandTree->GetBranch("muhat"))     bandTree->SetBranchAddress ("muhat",     &poihat1);
  if (bandTree->GetBranch("muhat_err")) bandTree->SetBranchAddress ("muhat_err", &poierr1);
  if (bandTree->GetBranch("tstype"))    bandTree->SetBranchAddress ("tstype",    &testStatType1);
  if (bandTree->GetBranch("wsfile"))    bandTree->SetBranchAddress ("wsfile",    &wsfile1);
  if (bandTree->GetBranch("optimize"))  bandTree->SetBranchAddress ("optimize",  &optimize1);
  if (bandTree->GetBranch("fitInRange"))bandTree->SetBranchAddress ("fitInRange",&fitInRangeBits);
  bandTree->GetEntry(0);
  if (wsfile1 && !wsfile1->empty()) {
    if       (wsfile.empty() || massAxis>=0)
      wsfile= *wsfile1;
    else if (*wsfile1!=wsfile) {
      cerr << "results from different workspace files: " << wsfile << " and " << *wsfile1 << endl;
      if (!force) return 2;
    }
  }
  if (invMass1>0.0) {
    if      (invMass<=0.0 || massAxis>=0)
      invMass= invMass1;
    else if (invMass1!=invMass) {
      cerr << "results at different masses: " << invMass << " and " << invMass1 << " GeV" << endl;
      if (!force) return 3;
    }
  }
  if (testStatType1>=0) {
    if (testStatType<0)
      testStatType= testStatType1;
    else if (testStatType1!=testStatType) {
      cerr << "results at different test statistic types: " << testStatType << " and " << testStatType1 << endl;
      if (!force) return 4;
    }
  }
  if (poierr==0.0 && poierr1!=0.0) {
    poihat= poihat1;
    poierr= poierr1;
  }
  optimize |= optimize1;

  if (fitInRangeBits != (UInt_t)-1) {
    if (fitInRange.empty()) {
      for (; fitInRangeBits; fitInRangeBits>>=1)
        fitInRange.push_back (fitInRangeBits&1);
      if (fitInRange.empty()) fitInRange.push_back(false);
    } else {
      UInt_t fitInRangeBits1=0;
      for (size_t i=0, n=fitInRange.size(); i<n; i++)
        if (fitInRange[i]) fitInRangeBits1 |= (1<<i);
      if (fitInRangeBits != fitInRangeBits1) {
        cerr << "results with different fit-in-range POIs: bit mask = " << fitInRangeBits1 << " and " << fitInRangeBits << endl;
        if (!force) return 5;
      }
    }
  }

  bandTree->ResetBranchAddresses();
  delete bandTree;
  return 0;
}


void
LikelihoodCalculator::PrintResult()
{
  cout << "Result";
  if (!wsfile.empty()) cout << " for workspace file " << wsfile;
  if (invMass>0.0) cout << " at " << invMass << " GeV";
  if (poierr!=0.0) cout << " with muhat = " << poihat << " +/- " << poierr;
  vector<Int_t> np= npoints;
  for (size_t i=0, n=np.size(); i<n; i++)
    if (fitInRange[i]) np[i]= -np[i];
  cout << " and " << Join(np,"x") << " " << Join (poiName) << " bins in range "<<Join(scanMin,scanMax);
  if (testStatType>=0) cout << " and test stat #" << testStatType;
  cout <<endl;
}


void LikelihoodCalculator::InterpolateMinimum (double xl, double x0, double xh, double yl, double y0, double yh,
                                               double &xm, double &ym)
{
  // Returns the minumum (xm,ym) of a quadratic interpolation of the
  // three points (xl,yl), (x0,y0), and (xh,yh), where xl<x0<xh.
  // We use the form: y(x) = (d2/2)*(x-x0)^2 + d1*(x-x0) + d0
  // Borrowed code from: http://www.ebyte.it/library/codesnippets/P3Interpolation.html
  double d2= 2.0 * ((yh-y0)/(xh-x0) - (yl-y0)/(xl-x0)) / (xh-xl);  // 2nd derivative
  double d1= (yh-y0)/(xh-x0) - 0.5*d2*(xh-x0);                     // 1st derivative at x0
  if (d2>0.0) {                  // minimum
    xm= x0 - d1/d2;
    ym= y0 + 0.5*d1*(xm-x0);
  } else if (yl==y0 && yh==y0) { // flat, so return (x0,y0)
    xm= x0; ym= y0;
  } else if (yh<yl)  {           // minimum at xh
    xm= xh; ym= yh;
  } else {                       // minimum at xl (or yl=yh)
    xm= xl; ym= yl;
  }
}


Double_t LikelihoodCalculator::FindMinimum (const TH1D* h, Int_t& imin, Double_t& xmin, int interpolate, int verbose)
{
  Double_t ymin= NaN;
  imin= -1;
  xmin= NaN;
  Int_t n=h->GetNbinsX();
  for (Int_t i=1; i<=n; i++) {
    Double_t v= h->GetBinContent(i);
    if (v!=0.0 && (imin==-1 || v<ymin)) {
      imin= i;
      ymin= v;
    }
  }
  vector<Double_t> points= GetBinPoints(h->GetXaxis());
  if (imin>0) xmin= points[imin-1];

  if (interpolate!=1 || imin<=1 || imin>=n) return ymin;
  if (h->GetBinContent(imin-1)==0.0)   return ymin;
  if (h->GetBinContent(imin+1)==0.0)   return ymin;

  double xint=xmin, yint=ymin;
  LikelihoodCalculator::InterpolateMinimum (points[imin-2],           xmin, points[imin],
                                            h->GetBinContent(imin-1), ymin, h->GetBinContent(imin+1),
                                            xint, yint);
  Int_t iint= h->FindFixBin(xint);
  if (iint != imin) {
    if (verbose>=0)
      cout << Form("Interpolation moved MLE from bin %d (%g,%g) to bin %d (%g,%g) at (%g,%g)",
                   imin, xmin, ymin, iint, points[iint-1], h->GetBinContent(iint), xint, yint) << endl;
  } else {
    if (verbose>=2)
      cout << Form("Interpolation moved MLE from bin %d (%g,%g) to (%g,%g)",
                   imin, xmin, ymin, xint, yint) << endl;
    ymin= yint;
    xmin= xint;
  }
  return ymin;
}


Int_t LikelihoodCalculator::FindErrors (const TH1D* h, Double_t xmin, const std::vector<Double_t>& up, TMatrixD& xval,
                                        int interpolate, int verbose)
{
  Int_t bin0= isnan(xmin) ? h->GetNbinsX()/2 : h->GetXaxis()->FindFixBin(xmin);
  int nup= up.size();
  xval.ResizeTo(0,-1,-nup,nup);  // initially 0 rows
  vector<Double_t> lxpm(2,NaN), lypm(2,-1.0);
  vector<Double_t> points= GetBinPoints(h->GetXaxis());
  for (Int_t istep=0, nb=h->GetNbinsX(); bin0-istep>=1 || bin0+istep<=nb; istep++)
    for (Int_t idir=0; idir<2; idir++) {           // look for UP points closest to minimum first (+,-,+,-,...)
      Int_t dir= 1-2*idir, i= bin0 + dir*istep;
      if (i<1 || i>nb) continue;
      Double_t x= points[i-1], y= h->GetBinContent(i);
      if (y==-1.0) continue;
      Double_t &lx=lxpm[idir], &ly=lypm[idir];   // remember last point for this direction
      if (ly!=-1.0) {
        for (Int_t j=0; j<nup; j++) {
          if (!((ly< up[j] && y>=up[j]) || (ly>=up[j] && y< up[j]))) continue;
          Double_t xup;
          if (interpolate==1) xup= lx + (x-lx) * (fabs(up[j])-ly) / (y-ly);
          else if (y<ly)      xup= lx;
          else                xup=  x;
          Int_t pm= dir * (y<ly ? -1 : 1);
          Int_t iup= pm*(j+1);
          Int_t imin=0, nmin=xval.GetNrows();
          for (; imin<nmin; imin++)
            if (isnan(xval(imin,iup))) break;
          if (imin>=nmin) {
            xval.ResizeTo (0,imin,-nup,nup);
            TMatrixDRow(xval,imin).Assign(NaN);
            xval(imin,0)= xmin;
          }
          xval(imin,iup)= xup;
        }
      }
      lx=x;
      ly=y;
    }

  Int_t nmin= xval.GetNrows();
  if (nmin>0 && verbose>=0) {
    cout << Form("%s = %.4f ", h->GetXaxis()->GetTitle(), xmin);
    for (int imin= 0; imin<nmin; imin++) {
      if (imin>0) cout << " and ";
      for (int j=0; j<nup; j++) {
        Double_t xlo= xval(imin,-j-1);
        Double_t xhi= xval(imin, j+1);
        if (isnan(xlo) && isnan(xhi)) continue;
        if (j>0) cout << " (";
        if      (isnan(xlo)) cout << Form("+%.4f [: %.4f]",xhi-xmin,xhi);
        else if (isnan(xhi)) cout << Form("-%.4f [%.4f :]",xmin-xlo,xlo);
        else {
          if (xlo<xmin && xhi>xmin) {
            Double_t errlo= xmin-xlo, errhi= xhi-xmin;
            if (fabs (errhi-errlo) < 0.01) cout << Form("+/- %.4f ",0.5*(errhi+errlo));
            else                           cout << Form("+%.4f -%.4f ",errhi,errlo);
          }
          cout << Form("[%.4f : %.4f]",xlo,xhi);
        }
        if (j>0) cout << ")";
      }
    }
    cout << endl;
  }

  return nmin;
}


LikelihoodCalculator::HistTypeBase*
LikelihoodCalculator::Profile (const HistTypeBase* h, Int_t axis, TObjArray* h1d,
                               TMatrixD* path, const TMatrixD* oldPath, const std::vector<Int_t>& axes,
                               int interpolate, int verbose)
{
  // Profiles over one variable (specified by "axis"), returning a new
  // histogram with one fewer dimensions.
  const TAxis* ax= h->GetAxis(axis);
  Int_t nd= h->GetNdimensions();
  vector<int> nbins(nd), nbprof(nd-1);
  vector<Double_t> xlo(nd-1), xhi(nd-1);
  int nscan=1;
  for (int d=0,dp=0; d<nd; d++) {
    nbins[d]= h->GetAxis(d)->GetNbins();
    if (d==axis) continue;
    nscan *= nbins[d];
    nbprof[dp]= nbins[d];
    xlo[dp]= h->GetAxis(d)->GetXmin();
    xhi[dp]= h->GetAxis(d)->GetXmax();
    dp++;
  }
  Int_t nax= ax->GetNbins();
  HistTypeBase* hp= new HistType (Form("%s_%s",h->GetName(),ax->GetTitle()),
                                  Form("%s profiled over %s",h->GetTitle(),ax->GetTitle()),
                                  nd-1, &nbprof[0], &xlo[0], &xhi[0]);
  for (int d=0,dp=0; d<nd; d++) {
    if (d==axis) continue;
    ReplaceAxis (hp->GetAxis(dp), h->GetAxis(d));
    dp++;
  }
  if (h1d) h1d->Expand(nscan);
  int ndpath= oldPath ? oldPath->GetNcols() : nd;
  vector<Int_t> pathAxes(ndpath,-1);
  if (oldPath) {
    for (int d=0, n=axes.size(); d<n; d++) pathAxes[axes[d]]= d;
  } else {
    for (int d=0; d<nd; d++)               pathAxes[d]= d;
  }
  if (path) path->ResizeTo(nscan,ndpath);
  for (int i=0; i<nscan; i++) {
    vector<int> binp= GetBin (i, nbprof, false, 1);
    TString binName= BinName(hp,binp);
    TH1D* hx= new TH1D (Form("%s_%s_%05d",h->GetName(),ax->GetTitle(),i),
                        Form("%s;%s;%s",binName.Data(),ax->GetTitle(),h->GetTitle()),
                        nax, ax->GetXmin(), ax->GetXmax());
    ReplaceAxis (hx->GetXaxis(), ax);
    hx->GetXaxis()->ResetAttAxis("X");  // THnD axes don't have attributes set
    hx->SetMinimum(0.0);
    hx->SetMarkerStyle(20);
    hx->SetMarkerSize(0.8);
    hx->SetMarkerColor(hx->GetLineColor());
    for (int j=1; j<=nax; j++) {
      vector<int> binh= GetBin (i, nbins, false, 1, axis, j);
      hx->SetBinContent (j, h->GetBinContent (&binh[0]));
    }
    Double_t xmin= NaN;
    Int_t    imin= -1;
    Double_t ymin= FindMinimum (hx, imin, xmin, interpolate, verbose);
    if (imin>0) {
      hp->SetBinContent (&binp[0], ymin);
      if (verbose>=1) {
        TH1D* hx2= dynamic_cast<TH1D*>(hx->Clone(Form("%s_off",hx->GetName())));
        hx2->SetTitle(Form(";At fixed %s, %s",binName.Data(),ax->GetTitle()));
        HistToLikelihood (hx2, ymin);
        TMatrixD xdummy;
        FindErrors (hx2, xmin, GetUpLevels(3,1), xdummy, interpolate, verbose);
        delete hx2;
      }
    }
    if (path) {
      vector<int> bin= GetBin (i,   nbins, false, 1, axis, imin);
      int        ibin= GetBin (bin, nbins, false, 1);
      for (int d=0; d<ndpath; d++) {
        int dd= pathAxes[d];
        (*path)(i,d)=  dd<0      ? (*oldPath)(ibin,d) :
                      (dd==axis) ? xmin
                                 : GetBinPoint(h,dd,bin[dd]);
      }
    }
    if (h1d) (*h1d)[i]= hx;
    else     delete hx;
  }
  return hp;
}


void LikelihoodCalculator::HistToLikelihood (TH1* h, Double_t mle)
{
  for (Int_t ix=1, nx=h->GetNbinsX(); ix<=nx; ix++)
    for (Int_t iy=1, ny=h->GetNbinsY(); iy<=ny; iy++)
      for (Int_t iz=1, nz=h->GetNbinsZ(); iz<=nz; iz++) {
        Int_t bin= h->GetBin(ix,iy,iz);
        Double_t v= h->GetBinContent(bin);
        Double_t nll= 2.0*(v-mle);
        h->SetBinContent (bin, v!=0.0 ? nll : -1.0);
      }
  if      (h->GetDimension()==1) h->GetYaxis()->SetTitle("#minus2ln#Lambda");
  else if (h->GetDimension()==2) h->GetZaxis()->SetTitle("#minus2ln#Lambda");
}


TH1D* LikelihoodCalculator::Profile2D (const HistTypeBase* h, int axis,
                                       TMatrixD* path, const TMatrixD* oldPath, const std::vector<Int_t>& axes,
                                       Double_t& mle, Double_t& xmle, Double_t& ymle, TMatrixD& xval)
{
  // For a 2D histogram, profiles over POI "axis", returning a 1D plot against the other POI.
  assert (h->GetNdimensions()==2);
  int axip= 1-axis;
  int interpolate1= interpolate && !fitInRange[axes[axis]];
  int interpolate2= interpolate && !fitInRange[axes[axip]];
  TObjArray hist1d;
  HistTypeBase* hpp= Profile (h, axis, &hist1d, path, oldPath, axes, interpolate1, verbose);
  TH1D* hp= hpp->Projection(0);
  hp->SetTitle(hpp->GetTitle());
  hp->GetXaxis()->SetLimits(hpp->GetAxis(0)->GetXmin(),hpp->GetAxis(0)->GetXmax());  // restore limits
  delete hpp;

  xmle= NaN;
  ymle= NaN;
  Double_t& xmle1= axis!=0 ? xmle : ymle;
  Double_t& xmle2= axis!=0 ? ymle : xmle;

  Int_t imle1= -1;
  mle= FindMinimum (hp, imle1, xmle1, interpolate2, verbose);
  if (imle1<=0) {
    hist1d.SetOwner();
    return hp;
  }

  xmle2= (*path)(imle1-1,axis);
  Int_t imle2= h->GetAxis(axis)->FindFixBin(xmle2);
  Int_t bin[2];
  bin[axip]= imle1;
  bin[axis]= imle2;
  cout << Form("%s,%s MLE is %.4f at %g,%g. Bin %d,%d is %.4f at %g,%g",
               h->GetAxis(axip)->GetTitle(), h->GetAxis(axis)->GetTitle(), mle, xmle1, xmle2,
               imle1, imle2, h->GetBinContent(bin),
               GetBinPoint(h,axip,imle1), GetBinPoint(h,axis,imle2)) << endl;
  HistToLikelihood (hp, mle);
  for (int i=0; i<hist1d.GetEntries(); i++)
    if (hist1d[i]) HistToLikelihood (dynamic_cast<TH1*>(hist1d[i]), mle);

  UseFittedPath (*path);

  FindErrors (hp, xmle1, upLevel, xval, interpolate, verbose);
  if (verbose>=0) PrintErrorPoints (hp, xval, path, upLevel);
  FixedErrors (hp, xval, hist1d);

  if (plotAllProjections) toPlot.AddAll (&hist1d);
  else                    hist1d.SetOwner();
  return hp;
}


void LikelihoodCalculator::FixedErrors (const TH1D* h, const TMatrixD& xval, const TObjArray& hist1d)
{
  for (int nupLevel=upLevel.size(), i=-nupLevel; i<=nupLevel; i++) {
    for (Int_t ival=0, nval=xval.GetNrows(); ival<nval; ival++) {
      Double_t xv=xval(ival,i);
      if (isnan(xv)) continue;
      const TAxis* ax= h->GetXaxis();
      Int_t  bin1= ax->FindFixBin(xv),   bin2=bin1;
      Double_t x1= ax->GetBinCenter(bin1), x2=x1;
      if (interpolate==1 && xv!=x1) {
        if (xv < x1)
          x1= ax->GetBinCenter(--bin1);
        else
          x2= ax->GetBinCenter(++bin2);
      }

      if (bin1<1 || bin2>hist1d.GetEntries()) continue;
      TH1D* hx1= dynamic_cast<TH1D*>(hist1d[bin1-1]);
      if (!hx1) continue;
      TH1D* hx2= dynamic_cast<TH1D*>(hist1d[bin2-1]);
      if (!hx2) continue;
      TH1D* hxi= dynamic_cast<TH1D*>(hx1->Clone(Form("%s_interp",hx1->GetName())));
      Double_t up= i ? upLevel[abs(i)-1] : 0.0;
      if (i<0) up= -up;
      hxi->SetTitle(Form("%s=%.4f (%s)",h->GetXaxis()->GetTitle(),xv,UpLevelName(up)));
      if (x2>x1) {
        Double_t xb= (xv-x1)/(x2-x1);
        for (Int_t j=1, n=hx1->GetNbinsX(); j<=n; j++) {
          Double_t y1= hx1->GetBinContent(j);
          Double_t y2= hx2->GetBinContent(j);
          hxi->SetBinContent (j, y1 + xb*(y2-y1));
        }
      }
      toPlot.Add(hxi);

      if (verbose>=0) {
        Double_t xmin= NaN;
        Int_t    imin= -1;
        Double_t ymin= FindMinimum (hxi, imin, xmin, interpolate, verbose);
        if (imin>0) {
          TH1D* hxo= dynamic_cast<TH1D*>(hxi->Clone(Form("%s_off",hxi->GetName())));
          hxo->SetTitle(Form(";At fixed %s, %s",hxi->GetTitle(),hxi->GetXaxis()->GetTitle()));
          for (Int_t j=1, n=hxi->GetNbinsX(); j<=n; j++)
            hxo->SetBinContent(j,hxi->GetBinContent(j)-ymin);
          TMatrixD xdummy;
          FindErrors (hxo, xmin, upLevel, xdummy, interpolate, verbose+1);
          delete hxo;
        }
      }
    }
  }
}


void LikelihoodCalculator::PrintErrorPoints (const TH1D* h, const TMatrixD& xval, TMatrixD* path, const std::vector<Double_t>& upLevel)
{
  cout << h->GetXaxis()->GetTitle() << " points at";
  int nu=0;
  for (int nupLevel=upLevel.size(), i=-nupLevel; i<=nupLevel; i++) {
    int ns=0;
    for (Int_t imin=0, nmin=xval.GetNrows(); imin<nmin; imin++) {
      Double_t xv=xval(imin,i);
      if (isnan(xv)) continue;
      if (ns==0) {
        if (nu>0) cout << ';';
        if (i==0) cout << " MLE=";
        else      cout << Form(" %+dsigma=",i);
      } else
        cout << " and ";
      Int_t bin= h->GetXaxis()->FindFixBin(xv);
      for (int d=0, n=path->GetNcols(); d<n; d++) {
        if (d>0) cout << ',';
        cout << Form("%.4f",(*path)(bin-1,d));
      }
      ns++;
      nu++;
    }
  }
  cout << endl;
}

void LikelihoodCalculator::PlotProfile (TPad* canvas, const TH1D* h, const TMatrixD& xval, Double_t mle, int axis)
{
  if (h->GetNbinsX() <= 1) return;
  TLine hline, vline;
  hline.SetLineWidth(1);
  hline.SetLineStyle(2);
  vline.SetLineWidth(2);

  TGraph g (h);
  vector<Double_t> points= GetBinPoints(h->GetXaxis());
  g.SetMarkerStyle(20);
  g.SetMarkerSize(0.8);
  g.SetLineColor(kBlack);
  g.SetLineWidth(2);
  g.SetMarkerColor(kBlack);
  // Remove missing points
  Double_t hmin= 0.0;  // keep minimum at 0
  for (int i=g.GetN()-1; i>=0; i--) {
    Double_t x=-1.0,y=-1.0;
    g.GetPoint(i,x,y);
    if (isnan(y) || y==-1.0)
      g.RemovePoint(i);
    else {
      g.SetPoint(i,points[i],y);
      if (y<hmin) hmin= y;
    }
  }
  TGraph* gcopy= dynamic_cast<TGraph*>(g.Clone(Form("graph_%s",poiName[axis].c_str())));
  gcopy->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  gcopy->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  saved.Add(gcopy);

  if (markerSM.size()==1) canvas->SetRightMargin(0.105);

  TH1D hc= *h;
  // Prevent TH1::Draw choking
  for (int i=1, n=hc.GetNbinsX(); i<=n; i++) {
    Double_t y= hc.GetBinContent(i);
    if (std::isinf(y) && y>0.0) hc.SetBinContent(i,1e120);
  }
  Double_t xlo=hc.GetXaxis()->GetXmin(), xhi=hc.GetXaxis()->GetXmax();
  hc.GetXaxis()->Set(hc.GetNbinsX(),xlo,xhi);  // override variable bins for axis plot
  hc.SetMinimum(hmin);
  if (hc.GetMaximum() > ymax) hc.SetMaximum(ymax);
  // Use a TGraph to draw points, so we don't show points >10 (TH1D plots them all at 10).
  hc.Draw("axis");
  g.Draw (showCaptions ? plotOpt : "c");
  static const Color_t color[]=  { kBlue, kRed, kGreen };
  for (int i=0, nupLevel=upLevel.size(); i<nupLevel; i++) {
    if (upLevel[i] > hc.GetMaximum()) break;
    if (i<int(sizeof(color)/sizeof(color[0]))) {
      hline.SetLineColor(color[i]);
      vline.SetLineColor(color[i]);
    }
    hline.DrawLine (xlo, upLevel[i], xhi, upLevel[i]);
    for (Int_t imin=0, nmin=xval.GetNrows(); imin<nmin; imin++) {
      for (int s=-1; s<=1; s+=2) {
        Double_t err=xval(imin,s*(i+1));
        if (isnan(err)) continue;
        vline.DrawLine (err, hmin, err, interpolate ? upLevel[i] : hc.GetBinContent(hc.FindFixBin(err)));
      }
    }
  }
  saved.Add (xval.Clone(), Form("limits_%s",poiName[axis].c_str()));

  TLatex t;
  t.SetTextColor(kBlack);
  t.SetTextAlign(12);
  t.SetTextSize(0.8*hc.GetYaxis()->GetLabelSize());
  t.SetTextFont(hc.GetYaxis()->GetLabelFont());
  double xs= xhi + hc.GetYaxis()->GetLabelOffset()*(xhi-xlo);
  int donesig= 0;
  if (markerSM.size()==1 && g.GetN()>0) {
    Double_t xsm= markerSM[0];
    Int_t pt;
    Double_t x1=xlo-1.0, y1=-1.0;
    for (pt=g.GetN()-1; pt>=0; pt--) {
      g.GetPoint(pt,x1,y1);
      if (x1<=xsm) break;
    }
    Double_t ysm;
    if (interpolate) {
      if      (pt<0)
        pt= 0;
      else if (pt>=g.GetN()-1) {
        pt= g.GetN()-2;
        g.GetPoint(pt,x1,y1);
      }
      Double_t x2, y2;
      g.GetPoint(pt+1,x2,y2);
      if (x2>x1) ysm= (xsm-x1)/(x2-x1)*(y2-y1) + y1;
      else       ysm= y1;
    } else
      ysm= y1;
    if (verbose>=0) {
      double Z= sqrt (ysm>0.0 ? ysm : 0.0);
      Double_t p= 2.0 * ROOT::Math::gaussian_cdf(-Z);
      cout << Form("At %s=%g: %s = %.4f", h->GetXaxis()->GetTitle(), xsm, h->GetYaxis()->GetTitle(), ysm)
           << (isnan(mle) ? "" : Form(" (%.4f)", mle+0.5*ysm))
           << Form(", p-value = %g, significance = %.2f (1DoF,2-sided)", p, Z)
           << endl;
    }
    if (xsm>=xlo && xsm<=xhi && ysm>hmin && ysm<=hc.GetMaximum()) {
      hline.SetLineColor(kOrange);
      vline.SetLineColor(kOrange);
      hline.DrawLine (xlo, ysm, xhi, ysm);
      vline.DrawLine (xsm, hmin, xsm, ysm);
      Double_t ylab= ysm;
      for (int i=0, nupLevel=upLevel.size(); i<nupLevel; i++) {
        donesig= i+1;
        Double_t dy= ylab-upLevel[i], yoff= 0.3-fabs(dy);
        if (yoff>0.0) {
          if (dy<0.0) yoff= -yoff;
          ylab= ysm+yoff;
          t.DrawLatex (xs, upLevel[i]-yoff, UpLevelName(upLevel[i]));
          break;
        } else
          t.DrawLatex (xs, upLevel[i],      UpLevelName(upLevel[i]));
      }
      t.DrawLatex (xs, ylab, Form("%s(%g)",h->GetYaxis()->GetTitle(),xsm));
    }
  }

  for (int i=donesig, nupLevel=upLevel.size(); i<nupLevel; i++)
    if (upLevel[i] <= hc.GetMaximum())
      t.DrawLatex (xs, upLevel[i], UpLevelName(upLevel[i]));

  TPaveText* label=0;
  if (showAtlasLabel) label= AtlasLabel (canvas, 0.35);

  PrintPlot (canvas);
  delete label;
  if (markerSM.size()==1) canvas->SetRightMargin(0.04);
}


TGraph* LikelihoodCalculator::Path2Graph (const TMatrixD& path, int axis, const std::vector<Int_t>& axes)
{
  int nm=0, nm0=path.GetNrows();
  TVectorD xm(nm0), ym(nm0);
  for (int i=0; i<nm0; i++) {
    if (isnan(path(i,axes[axis]))) continue;
    xm[nm]= path(i,axes[0]);
    ym[nm]= path(i,axes[1]);
    nm++;
  }
  if (!nm) return 0;
  xm.ResizeTo(nm);
  ym.ResizeTo(nm);
  TGraph* g= new TGraph (xm, ym);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->SetMarkerColor(g->GetLineColor());
  return g;
}


void LikelihoodCalculator::UseFittedPath (TMatrixD& path) const
{
  for (int d=0, nd=path.GetNcols(); d<nd; d++) {
    const HistTypeBase* h= res2[d];
    if (!h) continue;
    for (int i=0, n=path.GetNrows(); i<n; i++) {
      vector<Int_t> bin(nd);
      int dd;
      for (dd=0; dd<nd; dd++) {
        if (isnan(path(i,dd))) break;
        bin[dd]= h->GetAxis(dd)->FindFixBin(path(i,dd));
        if (bin[dd]<1 || bin[dd]>h->GetAxis(dd)->GetNbins()) break;
      }
      if (dd<nd) {
        if (verbose>=0 && !isnan(path(i,dd)))
          cerr << "MLE point at "<<Join(TMatrixDRow(path,i))<<" is outside histogram limits"<<endl;
        continue;
      }
      if (verbose>=2) cout << "Fitted path "<<Join(TMatrixDRow(path,i));
      path(i,d)= h->GetBinContent(&bin[0]);
      if (verbose>=2) cout << " -> "<<Join(TMatrixDRow(path,i))<<endl;
    }
  }
}


Double_t LikelihoodCalculator::Plot2D (TPad* canvas, const HistTypeBase* h, const std::vector<Int_t>& axesIn, TMatrixD* path)
{
  Int_t ndim= h->GetNdimensions();
  assert (ndim>=2);

  vector<Int_t> axes;
  if (axesIn.empty()) {
    axes.resize(ndim);
    for (Int_t d=0; d<ndim; d++) axes[d]= d;
  } else
    axes= axesIn;

  if (ndim>=3) {
    Double_t mle= NaN;
    for (Int_t d=0; d<ndim; d++) {
      TMatrixD profPath;
      HistTypeBase* hd= Profile (h, d, 0, &profPath, path, axes, interpolate && !fitInRange[axes[d]], verbose);
      vector<Int_t> newAxes;
      for (Int_t dd=0; dd<ndim; dd++)
        if (dd!=d) newAxes.push_back (axes[dd]);
      // Recursively profile each POI until we have 2D plots
      Double_t thisMle= Plot2D (canvas, hd, newAxes, &profPath);
      delete hd;
      if (isnan(mle)) mle= thisMle;
    }
    return mle;
  }

  Float_t oldoff= gStyle->GetTitleYOffset();
  if (atlasStyle) gStyle->SetTitleYOffset (atlasStyle>1 ? 1.6 : 1.75);
  TH1* h2= h->Projection(1,0);
  h2->SetTitle(h->GetTitle());
  gStyle->SetTitleYOffset (oldoff);
  //  h2->GetXaxis()->SetLimits(h->GetAxis(0)->GetXmin(),h->GetAxis(0)->GetXmax());  // restore limits
  //  h2->GetYaxis()->SetLimits(h->GetAxis(1)->GetXmin(),h->GetAxis(1)->GetXmax());  // restore limits
  if (atlasStyle>1 && h2->GetXaxis()->GetXmax()>=100.0)
    h2->GetXaxis()->SetNdivisions(50205);
  if (atlasStyle>1 && h2->GetYaxis()->GetXmax()>=100.0)
    h2->GetYaxis()->SetNdivisions(50205);

  TMatrixD path1, path2, xval1, xval2;
  Double_t mle1=NaN, mle2=NaN, xmle1=NaN, ymle1=NaN, xmle2=NaN, ymle2=NaN;
  TH1D* hp1= Profile2D (h, 1, &path1, path, axes, mle1, xmle1, ymle1, xval1);
  TH1D* hp2= Profile2D (h, 0, &path2, path, axes, mle2, xmle2, ymle2, xval2);
  if (!isnan(mle1)) HistToLikelihood (h2, mle1);

  TMarker pt (xmle1, ymle1, kMultiply);
  pt.SetMarkerSize(2);
  pt.SetMarkerColor(kBlue);
  TMarker sm (1.0, 1.0, kPlus);
  if (markerSM.size()==2) {
    sm.SetX(markerSM[0]);
    sm.SetY(markerSM[1]);
  } else {
    if (TString(h2->GetXaxis()->GetTitle()).BeginsWith("BR")) sm.SetX(0.0);
    if (TString(h2->GetYaxis()->GetTitle()).BeginsWith("BR")) sm.SetY(0.0);
  }
  sm.SetMarkerSize(2);
  sm.SetMarkerColor(kBlack);

  canvas->SetLeftMargin (atlasStyle>1 ? 0.166 : atlasStyle ? 0.14 : 0.12);
  if (h2->GetNbinsX() > 1 || h2->GetNbinsY() > 1) {
    canvas->SetRightMargin (atlasStyle ? 0.15 : 0.125);
    if (atlasStyle!=1) canvas->SetBottomMargin (atlasStyle ? 0.15 : 0.125);

    gStyle->SetPalette(1,0);
    h2->SetMaximum(1.05*upLevel2D.back());
    h2->SetContour (upLevel2D.size(), &upLevel2D[0]);
    h2->Draw("colz");

    pt.Draw();
    sm.Draw();

    TGraph* g1= Path2Graph (path1, 1, axes);
    if (g1) {
      g1->SetLineColor(kRed);
      g1->SetMarkerColor(kRed);
      g1->Draw(plotOpt);
    }
    TGraph* g2= Path2Graph (path2, 0, axes);
    if (g2) {
      g2->SetLineColor(kOrange);
      g2->SetMarkerColor(kOrange);
      g2->Draw(plotOpt);
    }

    PrintPlot (canvas);
    delete g2;
    delete g1;

    canvas->SetRightMargin(0.04);
    canvas->SetBottomMargin(gStyle->GetPadBottomMargin());
  }

  if (atlasStyle>1) canvas->SetBottomMargin(0.142);
  if (h2->GetNbinsX() > 2 && h2->GetNbinsY() > 2) {

    TString plotName= Form("%s_vs_%s",poiName[axes[1]].c_str(),poiName[axes[0]].c_str());
    TGraph* g0= new TGraph (1, &xmle1, &ymle1);
    g0->SetNameTitle (Form("contour_%s_0sigma",plotName.Data()), h2->GetTitle());
    pt.TAttMarker::Copy(*g0);
    g0->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
    g0->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
    saved.Add(g0);

    // Just plot 2 contours by default (change with -:3)
    Int_t styles[]= { 1, 2, 3, 4 };
    size_t ncont= min (size_t(ncontour2D), sizeof(styles)/sizeof(styles[0]));
    gStyle->SetPalette(ncont,styles);
    h2->SetContour (min (ncont, upLevel2D.size()), &upLevel2D[0]);
    h2->SetLineWidth(2);
    saved.Add(h2->Clone(Form("hist_%s",plotName.Data())));
    TH2* h2c= dynamic_cast<TH2*>(h2->Clone(Form("%s_contour",h2->GetName())));
    for (Int_t ix=1, nx=h2c->GetNbinsX(); ix<=nx; ix++)
      for (Int_t iy=1, ny=h2c->GetNbinsY(); iy<=ny; iy++)
        if (h2c->GetBinContent(ix,iy)==-1.0) h2c->SetBinContent(ix,iy,1e30);
    TCanvas ctmp("");   // temporary batch-mode canvas to receive contours
    ctmp.cd();
    h2c->Draw("cont list");
    gPad->Update();
    TList graphs;
    TObjArray graphs1(ncont);
    if (const TObjArray* contours= dynamic_cast<const TObjArray*>(gROOT->GetListOfSpecials()->FindObject("contours"))) {
      int nup=contours->GetSize();
      typedef vector<Double_t> vectorD;
      typedef vector<vectorD>  vectorD2;
      vectorD2 xlo(nup), xhi(nup), ylo(nup), yhi(nup);
      for (int iup=0; iup<nup; iup++) {
        int ng=0, np= 0;
        for (TIter it= dynamic_cast<const TList*>(contours->At(iup)); const TGraph* g= dynamic_cast<const TGraph*>(it());) {
          if (iup<int(ncont)) {
            TString gname= Form("contour_%s_%dsigma",plotName.Data(),iup+1);
            if (ng>0) gname += Form("_%d",ng);
            TGraph* gc= dynamic_cast<TGraph*>(g->Clone(gname));
            gc->SetLineStyle(styles[iup]);
            gc->SetLineColor(h2c->GetLineColor());
            gc->SetMarkerStyle(20);
            gc->SetMarkerSize(0.8);
            gc->SetMarkerColor(gc->GetLineColor());
            gc->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
            gc->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
            graphs.Add(gc);
            graphs1.AddAtAndExpand(gc,iup);
            saved.Add(gc);
            ng++;
          }
          if (verbose>=2) {
            for (int i=0, n=g->GetN(); i<n; i++) {
              vectorD xy(2,NaN);
              if (g->GetPoint(i,xy[0],xy[1])<0) continue;
              if (!np || xy[0]<xlo[iup][0]) xlo[iup]= xy;
              if (!np || xy[1]<ylo[iup][1]) ylo[iup]= xy;
              if (!np || xy[0]>xhi[iup][0]) xhi[iup]= xy;
              if (!np || xy[1]>yhi[iup][1]) yhi[iup]= xy;
              np++;
            }
          }
        }
      }
      if (verbose>=2) {
        cout << "TH2 countour " << h2c->GetXaxis()->GetTitle() << " extrema at";
        for (int iup=nup-1; iup>=0; iup--) cout << (xlo[iup].size()==2 ? Form(" %.4f,%.4f",xlo[iup][0],xlo[iup][1]) : " -");
        for (int iup=0;     iup<nup;iup++) cout << (xhi[iup].size()==2 ? Form(" %.4f,%.4f",xhi[iup][0],xhi[iup][1]) : " -");
        cout << endl;
        cout << "TH2 countour " << h2c->GetYaxis()->GetTitle() << " extrema at";
        for (int iup=nup-1; iup>=0; iup--) cout << (ylo[iup].size()==2 ? Form(" %.4f,%.4f",ylo[iup][0],ylo[iup][1]) : " -");
        for (int iup=0;     iup<nup;iup++) cout << (yhi[iup].size()==2 ? Form(" %.4f,%.4f",yhi[iup][0],yhi[iup][1]) : " -");
        cout << endl;
      }
    }
    canvas->cd();

    h2c->Draw("axis");
    for (TIter it= &graphs; TGraph* g= dynamic_cast<TGraph*>(it());) g->Draw("l");

    Double_t xlo=h2c->GetXaxis()->GetXmin(), xhi=h2c->GetXaxis()->GetXmax();
    Double_t ylo=h2c->GetYaxis()->GetXmin(), yhi=h2c->GetYaxis()->GetXmax();
    Double_t mlo=max(xlo,ylo), mhi=min(xhi,yhi);
    bool showsm= (sm.GetX()>=xlo && sm.GetX()<=xhi &&
                  sm.GetY()>=ylo && sm.GetY()<=yhi);
    bool m_vs_m= (fabs(0.5*(xlo+xhi)-125.0)<10.0 && fabs(0.5*(ylo+yhi)-125.0)<10.0);
    TLine dm0 (mlo,mlo,mhi,mhi);
    dm0.SetLineWidth(2);
    dm0.SetLineColor(kRed);
    pt.Draw();
    if (showsm) sm.Draw();
    if (m_vs_m) dm0.Draw();

    TPaveText* label=0;
    if (showAtlasLabel) label= AtlasLabel (canvas);

    Double_t xoff= 0.01+canvas->GetLeftMargin(), yoff= 0.99-canvas->GetTopMargin();
    TLegend leg (xoff+0.05, yoff-0.38, xoff+0.35, yoff-0.18);
    leg.SetBorderSize(0);
    leg.SetLineColor(0);
    leg.SetFillStyle(0);
    leg.SetFillColor(0);
    leg.AddEntry (&pt, "Best fit", "P");
    for (size_t i=0; i<ncont; i++) {
      TGraph* g= dynamic_cast<TGraph*>(graphs1.At(i));
      if (!g) continue;
      leg.AddEntry (g, Form("%s CL", UpLevelName(upLevel2D[i],2,2,0.0)), "L");
    }
    if (m_vs_m)   leg.AddEntry (&dm0, "#Deltam_{H}=0", "L");
    if (showsm) {
      if (m_vs_m) leg.AddEntry (&sm, "#Deltam_{H}=0 best fit", "P");
      else        leg.AddEntry (&sm, "SM",                     "P");
    }
    leg.Draw();

    PrintPlot (canvas);

    if (h2c->GetNbinsX()<30 && h2c->GetNbinsY()<70) {
      gStyle->SetPaintTextFormat(".2f");
      h2->Draw("text");
      for (TIter it= &graphs; TGraph* g= dynamic_cast<TGraph*>(it());) g->Draw("l");
      pt.Draw();
      if (showsm) sm.Draw();
      if (m_vs_m) dm0.Draw();
      PrintPlot (canvas);
      gStyle->SetPaintTextFormat("g");
    }

    gStyle->SetPalette(1,0);
    delete label;
    delete h2c;
  }
  canvas->SetLeftMargin(0.09);

  PlotProfile (canvas, hp1, xval1, mle1, axes[0]);
  PlotProfile (canvas, hp2, xval2, mle1, axes[1]);

  if (verbose>=1) FindExtrema (h, mle1, interpolate, verbose);

  delete hp1;
  delete hp2;
  delete h2;

  return mle1;
}


void
LikelihoodCalculator::DrawPlots (TPad* canvas, const TCollection& plots)
{
  if (plots.IsEmpty()) return;
  gStyle->SetOptTitle(1);
  canvas->SetTopMargin(0.07);
  for (TIter it= &plots; TObject* obj= it();) {
    TH1D* hx= dynamic_cast<TH1D*>(obj);
    if (!hx) continue;
    hx->Draw(plotOpt);
    PrintPlot (canvas);
  }
  gStyle->SetOptTitle (0);
  canvas->SetTopMargin(0.04);
}

void LikelihoodCalculator::FindExtrema (const HistTypeBase* h, Double_t mle, int interpolate, int verbose, int nup)
{
  Double_t badmle, mlefac;
  if (mle==0.0) { // already converted to PLR
    badmle= -1.0; mlefac= 1.0;
  } else {
    badmle=  0.0; mlefac= 2.0;
  }
  int ndim= h->GetNdimensions();
  const vector<Double_t> upLevelND= GetUpLevels (nup, ndim, verbose);
  typedef vector<int>      vectorI;
  typedef vector<vectorI>  vectorI2;
  typedef vector<vectorI2> vectorI3;
  vectorI3 lo (ndim, vectorI2 (nup, vectorI(ndim)));
  vectorI3 hi (ndim, vectorI2 (nup, vectorI(ndim)));
  for (int i=0, n=h->GetNbins(); i<n; i++) {
    vector<Int_t> bin(ndim);
    Double_t v= h->GetBinContent(i,&bin[0]);
    if (mle==badmle) continue;
    v= mlefac*(v-mle);
    for (int j=0; j<nup; j++) {
      if (v>upLevelND[j]) continue;
      for (int d=0; d<ndim; d++) {
        if (!lo[d][j][d] || bin[d]<lo[d][j][d]) lo[d][j]= bin;
        if (!hi[d][j][d] || bin[d]>hi[d][j][d]) hi[d][j]= bin;
      }
    }
  }

  for (int d=0; d<ndim; d++) {
    TAxis* ax= h->GetAxis(d);
    cout << Form("%-10s",ax->GetTitle()) << " "<<-nup<<" to "<<nup<<"sigma at ";
    for (int dd=0; dd<ndim; dd++) {
      if (dd>0) cout << ',';
      cout << h->GetAxis(dd)->GetTitle();
    }
    cout << " =";
    int nb= ax->GetNbins();
    for (int s=-nup; s<=nup; s++) {
      if (!s) continue;
      int j= abs(s)-1;
      vectorI& bin= s<0 ? lo[d][j] : hi[d][j];
      TString out;
      if (bin[d] > 1 && bin[d] < nb) {  // can't tell at edges
        if (verbose>=1) out= Join(bin)+" = ";
        for (int dd=0; dd<ndim; dd++) {
          if (dd>0) out += ',';
          Double_t x= GetBinPoint(h,dd,bin[dd]);
          if (interpolate==1 && dd==d) {
            vector<Int_t> bin2= bin;
            for (int sgn=s<0?-1:1, k=bin[d]+sgn; k>=1 && k<=nb; k+=sgn) {
              bin2[d]= k;
              Double_t y2= h->GetBinContent(&bin2[0]);
              if (y2 == badmle) continue;
              y2= mlefac*(y2-mle);
              Double_t y= mlefac*(h->GetBinContent(&bin[0])-mle), x2= GetBinPoint(ax,k);
              Double_t xi= x2 + (x-x2) * (upLevelND[j]-y2) / (y-y2);
              if (verbose>=1 && fabs(xi-x)>0.01) out += Form("(%.4f,%.4f,%.4f,%.4f)->",x,y,x2,y2);
              x= xi;
              break;
            }
          }
          out += Form("%.4f",x);
        }
      }
      else
        out= "-";
      if (verbose>=1) cout << "  "<<out;
      else            cout << Form(" %-14s",out.Data());
    }
    cout << endl;
  }
}


void
LikelihoodCalculator::PlotResult(TPad* canvas)
{
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.04);

  res->SetTitle("#minus2ln#Lambda");   // all plots will be of -ln(PLR), so set the title here so they'll inherit it

  Double_t mle= NaN;
  if (res->GetNdimensions()>=2) {
    mle= Plot2D (canvas, res);
  } else {
    TH1D* h1= dynamic_cast<TH1D*>(res->Projection(0));
    h1->SetTitle(res->GetTitle());
    h1->GetXaxis()->SetLimits(res->GetAxis(0)->GetXmin(),res->GetAxis(0)->GetXmax());  // restore limits
    Double_t xmle= NaN;
    Int_t imle= -1;
    mle= FindMinimum (h1, imle, xmle, interpolate && !fitInRange[0], verbose);
    if (imle>0) {
      cout << Form("%s MLE is %.4f at %g. Bin %d is %.4f at %g",
                   h1->GetXaxis()->GetTitle(), mle, xmle,
                   imle, h1->GetBinContent(imle), GetBinPoint(h1->GetXaxis(),imle)) << endl;
      HistToLikelihood (h1, mle);
    }
    TMatrixD xval;
    FindErrors (h1, xmle, upLevel, xval, interpolate, verbose);
    PlotProfile (canvas, h1, xval, mle);
    delete h1;
  }

  if (verbose>=1) FindExtrema (res, mle, interpolate, verbose);

  // All saved plots at the end
  DrawPlots(canvas,toPlot);

  for (int i=0; i<nscan; i++) {
    vector<int> bin= GetBin (i, npoints, true, 1);
    Double_t v= res->GetBinContent(&bin[0]);
    if (v==0.0) {
      if (verbose>=0)
        cout << "Point "<<i<<" ("<<Join(bin)<<") at "<<BinName(res,bin)<< " is missing" << endl;
    } else if (verbose>=2)
      cout << "Point "<<i<<" at "<<BinName(res,bin)<< " is " << 2.0*(v-mle) << endl;
  }
}




//==============================================================================
// Constructors and destructor
//==============================================================================

LikelihoodCalculator::LikelihoodCalculator (const char* name)
  : WorkspaceCalculator (name)
{
  Init();
}

LikelihoodCalculator::LikelihoodCalculator (const char* name, int argc, const char* const* argv)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (argc, argv);
  cmdArgs= Join(vector<string>(argv+1,argv+argc)," ");
}

LikelihoodCalculator::LikelihoodCalculator (const char* name, const char* args)
  : WorkspaceCalculator (name)
{
  Init();
  initErr= ParseArgs (args);
  cmdArgs= args;
}

LikelihoodCalculator::~LikelihoodCalculator()
{
}


void LikelihoodCalculator::Usage() const
{
  cout << "Usage:"
       << "\n  "<< GetName()<<" [ WORKSPACE-FILE.root | RESULT-FILES.root ] \\"
       << "\n      -w workspace|'combWS' -m modelConfig|'"<<modelSBName<<"' -d dataset|'"<<dataName<<"' -S Seed|FILE.root[:N] \\"
       << "\n      -a (all points) -p jobNum -n nPointsToScan|"<<Join(npoints,"x")<<" -N points/job -P pointsToScan \\"
       << "\n      -M invMass -z scanMin:scanMax|"<<Join(scanMin,scanMax)<<" -2 poimax -c poiName|"<<Join(poiName,":")<<" \\"
       << "\n      -O opt-level|";
  PrintOptimize(optimize,0);
  cout
       <<" -T testStatType|"<<testStatType<<" -e fitTol -H (plot) \\"
       << "\n      -v (verbose) -q (quiet) -r OUTPUT-RESULT-FILE"
       << endl;
}

int LikelihoodCalculator::ParseArgs (int argc, const char* const* argv)
{
  const char*  pointsToScanArg = 0;
  bool         all             = false;
  int          nPointsNow      = -1;
  vector<string> optimizeStrings;

  for (int i= 1; i<argc; i++) {
    if (argv[i][0]!='-') {
      fileName.push_back(argv[i]);
      continue;
    }
    for (const char* a= argv[i]+1; *a;) {
      switch (const char c= *a++) {
        case 'h': case '?':                    break;
        case 'H': plotResult++;                break;
        case 'v': verbose++;                   break;
        case 'q': verbose--;                   break;
        case 'D': detailedOutput++;            break;
        case 't': doStatsTree--;               break;
        case 'I': interpolate=0;               break;
        case 'A': plotAllProjections=1;        break;
        case 'f': force = true;                break;
        case '.': plotOpt="c";                 break;
        default : if (!*a && i+1<argc) a= argv[++i];
        switch (c) {
        case 'w': wsName          = a;
                  altResultName   = true;      break;
        case 'm': modelSBName     =        a;  break;
        case 'd': dataName        =        a;  break;
        case 'r': resultFileName  =        a;  break;
        case 'P': pointsToScanArg =        a;  break;
        case 'S': seedName        =        a;  break;
        case 'L': initialSnapshot =        a;
                  optimize |= kLoadInitialSnapshot;
                  break;
        case 'O': optimizeStrings.push_back(a);break;
        case 'C': a= Scan (a, poiName, ':');   break;
        case 'x': a= Scan (a, extraPOI);       break;
#ifdef USE_ProfileLikelihoodTestStatEnhanced
        case 'l': initialSnapshotTS=       a;  break;
        case 'o': if (optimizeTS.Length()) optimizeTS += "+";
                  optimizeTS     +=        a;  break;
#endif
        default : const char* ai= a;
        switch (c) {
        case 'p': jobNum          = Strtol(a); break;
        case 'n': a= Scan (a, npoints, 'x', 1);break;
        case 'N': nPointsNow      = Strtol(a);
                  all             = false;     break;
        case 'M': invMass         = Strtod(a); break;
        case 'u': a= Scan (a, bPOI, ':');
                  if (*a==',') a= Scan (a+1, sbPOI, ':');
                  break;
        case '0': a= Scan (a, scanMin);        break;  // obsolete
        case '1': a= Scan (a, scanMax);        break;  // obsolete
        case '2': a= Scan (a, poimax);         break;
        case 'z': a= Scan (a, scanMin, scanMax);
                  poimin= scanMin;
                  if (poimax.empty()) poimax= scanMax;
                  break;
        case 'Z': a= Scan (a, *varSettings);   break;
        case 'T': testStatType    = Strtol(a); break;
        case 'e': if (*a=='p' && *++a) {
                    fitPrec= Strtod(a);
                    if (*a==',') a++;
                  }
                  if (*a && *a!=',') fitTol= Strtod(a); // Target EDM=0.001*fitTol. Default RooFit fitTol=1 (overrides Minuit fitTol=0.01).
                  if (*a==',' && *++a) {
                    fitTol0= fitTol;
                    fitTol= Strtod(a);
                  }
                  break;
        case 's': skipPlot        = Strtol(a); break;
        case '~': asimovSig       = Strtod(a); break;
        case 'F': maxFunctionCalls= Strtol(a); break;
        case 'a': newParameterRanges= Strtod(a);
                  optimize |= kAdjustRanges;   break;
        case 'E': fixStatErrorMin = Strtod(a); 
                  optimize |= kFixStatError;   break;
        case 'i': if (*a==',')         interpCodeNorm= -1;
                  else                 interpCodeNorm=  Strtol(a);
                  if (*a==',' && *++a) interpCodeShape= Strtol(a);
                  else                 interpCodeShape= -1;
                  optimize |= kInterpolateErrors;
                  break;
        case 'X': a= Scan (a, markerSM, ':');  break;
        case 'y': ymax= Strtod(a);             break;
        case ':': nupLevels= Strtol(a);
                  if (nupLevels<0) { lim95= true; nupLevels= -nupLevels; }
                  ncontour2D= nupLevels;
                  break;
        default:          cerr << argv[0] << ": invalid option -"       << c              << endl; return 1; }
        if      (*a)    { cerr << argv[0] << ": invalid option value -" << c << " " << ai << endl; return 1; }
        else if (a==ai) { cerr << argv[0] << ": missing -" << c << " option value"        << endl; return 1; }}
        a= ""; // go to next arg
      }
    }
  }

  for (size_t i= 0, n=optimizeStrings.size(); i<n; i++)
    if (!SetOptimize(optimizeStrings[i].c_str(),"-O")) return 1;

  if (seedName.Length()>0) {
    char* a= 0;
    long int setSeed = strtol (seedName.Data(), (char**)&a, 10);
    if (*a == '\0') {
      seed= setSeed;
      seedName= "";
    }
  }

  if (pointsToScanArg) {  // -P 1,2,3,4,10/2,/10 -> 1,2,3,4,6,8,10,20,30,...,scanMax
    vector<string> args;
    Scan (pointsToScanArg, args, ':');
    npoi= args.size();
    if (int(npoints.size()) > npoi) npoi= npoints.size();
    else                            npoints.resize(npoi,1);
    pointsToScan.resize(npoi);
    args.resize(npoi);
    for (size_t d=0; d<size_t(npoi); d++) {
      vector<double> v, s;
      if (*Scan (args[d].c_str(), v, s, '/', ',', true)) {
        cerr << argv[0] << ": invalid option value -P " << pointsToScanArg << endl;
        return 1;
      }
      if (v.empty()) continue;
      Double_t scanMin0= d<scanMin.size() ? scanMin[d] : NaN;
      Double_t scanMax0= d<scanMax.size() ? scanMax[d] : NaN;
      Double_t lastp= scanMin0;
      for (size_t i=0, n=v.size(); i<n; i++) {
        Double_t pointVal= v[i];
        if (!(isnan(s[i]) && std::signbit(s[i]))) {  // has / (not a single number)
          if (isnan(pointVal)) pointVal= scanMax0;
          Double_t step= s[i];
          if (!(step > 0.0) || isnan(pointVal) || isnan(lastp)) {  // step<=0 or nan
            cerr << argv[0] << ": invalid step size in -P " << pointsToScanArg << endl;
            return 1;
          }
          Double_t hi= pointVal + 1e-6*step;
          Double_t p= lastp + step;
          for (; p < hi; p += step)
            if (!(p<scanMin0 || p>=scanMax0+1e-6*step))    // scanMin/Max0=NaN if no limit
              pointsToScan[d].push_back(p>=scanMax0?scanMax0:p);
        } else if (!(pointVal<scanMin0 || pointVal>scanMax0))
          pointsToScan[d].push_back(pointVal);
        lastp= pointVal;
      }
      npoints[d]= pointsToScan[d].size();
    }
  } else if (!npoints.empty()) {
    npoi= npoints.size();
    pointsToScan.resize(npoi);
  }

  if (!npoints.empty()) {
    // -n n1,-n2 fits mu2 in the range of each bin.
    fitInRange.resize(npoints.size());
    for (size_t i=0, n=npoints.size(); i<n; i++) {
      if ((fitInRange[i]= (npoints[i]<0))) npoints[i]= -npoints[i];
    }

    nscan= 1;
    for (int i=0; i<npoi; i++) nscan *= npoints[i];
    if (nscan<=0) {
      cerr << argv[0] << ": bad number of points: " << Join(npoints,"x") << endl;
      return 1;
    }
  }

  if (jobNum<0 && nPointsNow<=0) all = true;
  if (jobNum<0) jobNum = 0;
  jobSet= jobNum;
  if (!all) {
    if (nPointsNow<=1) {
      firstPoint = lastPoint = jobNum % nscan;
      jobSet                 = jobNum / nscan;
    } else {
      jobNum *= nPointsNow;
      int np = ((nscan+nPointsNow-1)/nPointsNow)*nPointsNow;
      firstPoint = jobNum % np;
      jobSet     = jobNum / np;
      lastPoint = firstPoint+nPointsNow-1;
      if (lastPoint>=nscan) lastPoint = nscan-1;
    }
  }

  upLevel=   GetUpLevels (nupLevels, 1, 1, lim95);
  upLevel2D= GetUpLevels (6, 2);   // ncontour2D is just for plotting the contours
  saved.Add (new TVectorD(nupLevels,&upLevel[0]), "upLevels");

  if (fileName.empty()) {
    SetDefaults();
    Usage();
    return 1;
  }

  return 0;
}
