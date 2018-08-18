// $Id: limitPlot.C 503374 2012-05-31 18:27:22Z adye $
// Plot Higgs cross-section upper limits at different masses.
// Uses results from StandardHypoTestInv.C.
#include <iostream>
#include <iomanip>

#include "TString.h"
#include "TPRegexp.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TObject.h"
#include "TList.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"

#include "AtlasStyle.C"
#include "AtlasLabels.C"

using std::cout;
using std::endl;
using std::setw;

TCanvas* c1=0;

void printLimit (double lim, double err=0, double muMin=0, double muMax=1e30, bool showerr=true)
{
  bool errok= err>0.0 && err<1e30;
  cout << std::fixed << setprecision(3) << setw(7) << lim;
  if (showerr &&  errok) cout << " +/- " << setw(4) << err;
  cout << (muMax>muMin && (lim<=muMin||lim>=muMax)?'*':' ');
  if (showerr && !errok) cout << setw(10) << " ";
}

bool SetBranchAddress (TTree& t, const char* k, Double_t* v) { return t.GetBranch(k) ? (t.SetBranchAddress(k,v),true) : false; }
bool SetBranchAddress (TTree& t, const char* k,    Int_t* v) { return t.GetBranch(k) ? (t.SetBranchAddress(k,v),true) : false; }

void limitPlot (int         logmu= 1,   // 0=linear, 1=log plots, 2=suppress error printout, 4=log M
                const char* files= "Results*.root",  // specify multiple files with wildcards and/or separated by spaces
                double      mlo=   118.0,
                double      mhi=   192.0,
                double      mulo=  0,
                double      muhi= -1,
                int         pcl=   -99,   // 0 or -1 if !useCLs
                const char* compare= 0,
                const char* compare2= 0,
                const char* cmpname= 0,
                const char* cmp2name= 0)
{
  bool showerr= !(logmu&2), showerr_obs=showerr;
  if (muhi<=mulo) {
    if (logmu&1) {
      mulo=  0.1;
      muhi=  1e3;
    } else {
      mulo=  0;
      muhi=  40.5;
    }
  }

  SetAtlasStyle();

  cout << "Read limit files " << files << endl;
  TChain bandTree ("band", "limitPlot");
  TStringToken filespec (files, "\\s+");
  while (filespec.NextToken()) bandTree.Add (filespec);
  Int_t nm = bandTree.GetEntries();
  TIter next(bandTree.GetListOfFiles());
  const TChainElement* ntfile;
  TString nam;
  bool samenam= true;
  TPMERegexp re("(?:^|/)H_([^_*?/]+)_[1-6][0-9][05][_.-]");
  while ((ntfile= dynamic_cast<const TChainElement*>(next()))) {
    cout << "Limit file " << ntfile->GetTitle() << " has " << ntfile->GetEntries() << " entries" << endl;
    if  (samenam && re.Match(ntfile->GetTitle()) && (nam.Length()==0 || re[1]==nam)) {
      nam= re[1];
    } else {
      samenam= false;
      nam= "";
    }
  }
  if (nam.Length()>0) cout << "Samples " << nam << endl;

  double invMass=0;
  double band2sigDown=0, band1sigDown=0, bandMedian=0, band1sigUp=0, band2sigUp=0;
  double band2sigDown_err=0, band1sigDown_err=0, bandMedian_err=0, band1sigUp_err=0, band2sigUp_err=0;
  double obsUL=0, obsUL_err=0, CLb=0, muMin=0, muMax=0, confidenceLevel=0;
  Int_t nPointsToScan=0, nToysSig=0, nToysBkg=0;
  Bool_t useCLs=(pcl<-2), haveCLs=false, haveCLb=false;
  bandTree.SetBranchAddress ("invMass",      &invMass);
  bandTree.SetBranchAddress ("band2sigDown", &band2sigDown);
  bandTree.SetBranchAddress ("band1sigDown", &band1sigDown);
  bandTree.SetBranchAddress ("bandMedian",   &bandMedian);
  bandTree.SetBranchAddress ("band1sigUp",   &band1sigUp);
  bandTree.SetBranchAddress ("band2sigUp",   &band2sigUp);
  bandTree.SetBranchAddress ("obsUL",        &obsUL);
  bandTree.SetBranchAddress ("muMin",        &muMin);
  bandTree.SetBranchAddress ("muMax",        &muMax);
  bandTree.SetBranchAddress ("confidenceLevel",&confidenceLevel);
  bandTree.SetBranchAddress ("nPointsToScan",&nPointsToScan);
  bandTree.SetBranchAddress ("nToysSig",     &nToysSig);
  bandTree.SetBranchAddress ("nToysBkg",     &nToysBkg);
  if (bandTree.GetBranch("obsUL_err"))
    bandTree.SetBranchAddress ("obsUL_err",        &obsUL_err);
  else
    showerr_obs= false;
  if (bandTree.GetBranch("CLb")) {
    haveCLb= true;
    bandTree.SetBranchAddress ("CLb",          &CLb);
  }
  if (bandTree.GetBranch("bandMedian_err")) {
    bandTree.SetBranchAddress ("band2sigDown_err", &band2sigDown_err);
    bandTree.SetBranchAddress ("band1sigDown_err", &band1sigDown_err);
    bandTree.SetBranchAddress ("bandMedian_err",   &bandMedian_err);
    bandTree.SetBranchAddress ("band1sigUp_err",   &band1sigUp_err);
    bandTree.SetBranchAddress ("band2sigUp_err",   &band2sigUp_err);
  } else
    showerr= false;
  if (bandTree.GetBranch("useCLs")) {
    bandTree.SetBranchAddress ("useCLs",           &useCLs);
    haveCLs= true;
  }
    

  int w= 9+(showerr?9:0);
  int v= 9+(showerr_obs?9:0);
  cout << setw(9) << "M (GeV)"
       << setw(w) << "-2sigma"
       << setw(w) << "-1sigma"
       << setw(w) << "exp"
       << setw(w) << "+1sigma"
       << setw(w) << "+2sigma"
       << setw(v) << "obs";
  if (haveCLb)
  cout << setw(9) << "CLb";
  cout << setw(9) << "mu Step"
       << setw(9) << "mu Max"
       << setw(9) << "#toy S+B"
       << setw(9) << "#toy B";
  cout << endl;
  TVectorD x(nm), c(nm), obsCLsb(nm), obs(nm), ebp1(nm), ebm1(nm), ebp2(nm), ebm2(nm), zero(nm);
  Double_t* lim[]= {&band2sigDown, &band1sigDown, &bandMedian, &band1sigUp, &band2sigUp};
  size_t nlim= sizeof(lim)/sizeof(lim[0]);
  for (int i=0;i<nm;i++){
    bandTree.GetEntry(i);
    Double_t step= (nPointsToScan>0 ? (muMax-muMin)/nPointsToScan : 0);
    cout << std::fixed << setprecision(0)
         << setw(8) << invMass << ' ';
    printLimit (band2sigDown, band2sigDown_err, muMin, muMax, showerr);
    printLimit (band1sigDown, band1sigDown_err, muMin, muMax, showerr);
    printLimit (bandMedian,   bandMedian_err,   muMin, muMax, showerr);
    printLimit (band1sigUp,   band1sigUp_err,   muMin, muMax, showerr);
    printLimit (band2sigUp,   band2sigUp_err,   muMin, muMax, showerr);
    printLimit (obsUL,        obsUL_err,        muMin, muMax, showerr_obs);
    if (haveCLb)
    cout << std::fixed << setprecision(2)
         << setw(8) << CLb << ' ';
    cout << std::fixed << setprecision(2)
         << setw(8) << step << ' '
         << setw(8) << muMax << ' '
         << setw(8) << nToysSig << ' '
         << setw(8) << nToysBkg;
    if (haveCLs)
    cout << (useCLs ? " CLs" : " PCL");
    cout << endl;
    double observed= obsUL;
    if (!useCLs && pcl>=-2 && pcl<=2 && nlim==5) {
      if (observed < *lim[2+pcl]) observed= *lim[2+pcl];
    }
    if (muMax>muMin) {
      for (size_t j=0; j<nlim; j++) {
        if (*lim[j]<muMin) { *lim[j]= mulo; }
        if (*lim[j]>muMax) { *lim[j]= muhi; }
      }
    }
    x[i]= invMass;
    obsCLsb[i]= obsUL;
    obs[i]= observed;
    c[i]= bandMedian;
    ebp1[i]=band1sigUp-bandMedian;
    ebm1[i]=bandMedian-band1sigDown;
    ebp2[i]=band2sigUp-bandMedian;
    if (useCLs)
      ebm2[i]=bandMedian-band2sigDown;
    else
      ebm2[i]=bandMedian-band1sigDown;  // -1sigma from PCL
    zero[i]= 0.0;
  }
  bandTree.ResetBranchAddresses();


  TVectorD cmpx, cmpc, cmpm2s, cmpm1s, cmpp2s, cmpp1s, cmpobs;
  bool haveCmpBands= 0;
  if (compare) {
    double invMassCmp=0, band2sigDownCmp=0, band1sigDownCmp=0, bandMedianCmp=0, band1sigUpCmp=0, band2sigUpCmp=0, obsULCmp=0;
    cout << "Read comparison files " << compare << endl;
    TChain cmpTree ("band", "limitPlot");
    cmpTree.Add (compare);
    Int_t nc = cmpTree.GetEntries();
    cmpx.ResizeTo(nc);
    cmpc.ResizeTo(nc);
    cmpobs.ResizeTo(nc);
    cmpm2s.ResizeTo(nc);
    cmpm1s.ResizeTo(nc);
    cmpp2s.ResizeTo(nc);
    cmpp1s.ResizeTo(nc);
    cmpTree.SetBranchAddress ("invMass",      &invMassCmp);
    cmpTree.SetBranchAddress ("band2sigDown", &band2sigDownCmp);
    cmpTree.SetBranchAddress ("band1sigDown", &band1sigDownCmp);
    cmpTree.SetBranchAddress ("bandMedian",   &bandMedianCmp);
    cmpTree.SetBranchAddress ("band1sigUp",   &band1sigUpCmp);
    cmpTree.SetBranchAddress ("band2sigUp",   &band2sigUpCmp);
    cmpTree.SetBranchAddress ("obsUL",        &obsULCmp);
    for (int i=0;i<nc;i++){
      cmpTree.GetEntry(i);
      cmpx[i]= invMassCmp;
      cmpobs[i]= obsULCmp;
      cmpc[i]= bandMedianCmp;
      cmpp2s[i]= band2sigUpCmp;
      cmpp1s[i]= band1sigUpCmp;
      cmpm2s[i]= band2sigDownCmp;
      cmpm1s[i]= band1sigDownCmp;
      if (band2sigUpCmp || band2sigDownCmp) haveCmpBands= 2;
      else if               (bandMedianCmp) haveCmpBands= 1;
    }
  }

  TVectorD cmp2x, cmp2c, cmp2obs, cmp2err;
  if (compare2) {
    double invMassCmp=0, bandMedianCmp=0, obsULCmp=0, obsULCmp_err=0;
    cout << "Read 2nd comparison files " << compare2 << endl;
    TChain cmpTree ("band", "limitPlot");
    cmpTree.Add (compare2);
    Int_t nc = cmpTree.GetEntries();
    cmp2x.ResizeTo(nc);
    cmp2c.ResizeTo(nc);
    cmp2obs.ResizeTo(nc);
    cmpTree.SetBranchAddress ("invMass",      &invMassCmp);
    cmpTree.SetBranchAddress ("bandMedian",   &bandMedianCmp);
    cmpTree.SetBranchAddress ("obsUL",        &obsULCmp);
    if (SetBranchAddress (cmpTree, "obsUL_err",  &obsULCmp_err))  cmp2err.ResizeTo(nc);
    for (int i=0;i<nc;i++){
      cmpTree.GetEntry(i);
      cmp2x[i]= invMassCmp;
      cmp2obs[i]= obsULCmp;
      cmp2c[i]= bandMedianCmp;
      if (cmp2err.GetNrows()>0) cmp2err[i]= obsULCmp_err;
    }
  }

  TGraph* go= 0;
  if (!useCLs) {
    go= new TGraph (x, obsCLsb);
    go->SetLineColor(1);
    go->SetLineWidth(2);
    go->SetLineStyle(3);
  }

  TGraph* gp= new TGraph (x, obs);
  gp->SetLineColor(1);
  gp->SetLineWidth(2);
  gp->SetMarkerStyle(20);
  gp->SetMarkerColor(1);

  TGraph* gc= new TGraph (x, c);
  gc->SetLineColor(kBlack);
  gc->SetLineStyle(2);
  gc->SetLineWidth(2);

  TGraph* gb1= new TGraphAsymmErrors (x, c, zero, zero, ebm1, ebp1);
  gb1->SetFillColor(3);

  TGraph* gb2= new TGraphAsymmErrors (x, c, zero, zero, ebm2, ebp2);
  gb2->SetFillColor(5);

  TGraph *gco=0, *gcc=0, *gcp1s=0, *gcm1s=0;
  if (compare) {
    gco= new TGraph (cmpx, cmpobs);
    gco->SetLineWidth(2);
    gco->SetLineColor(4);
    gco->SetMarkerColor(kRed);
    //    gco->SetMarkerSize(1.2);
    gco->SetMarkerStyle(24);

    if (haveCmpBands>=1) {
      gcc= new TGraph (cmpx, cmpc);
      gcc->SetLineStyle(2);
      gcc->SetLineWidth(2);
      gcc->SetLineColor(4);
    }

    if (haveCmpBands>=2) {
      gcp1s= new TGraph (cmpx, cmpp1s);
      gcp1s->SetLineStyle(2);
      gcp1s->SetLineWidth(2);
      gcp1s->SetLineColor(4);

      gcm1s= new TGraph (cmpx, cmpm1s);
      gcm1s->SetLineStyle(2);
      gcm1s->SetLineWidth(2);
      gcm1s->SetLineColor(4);
    }
  }

  TGraph *gc2o=0, *gc2c=0;
  if (compare2) {
    if (cmp2err.GetNrows()>0) {
      TVectorD zero2(cmp2x.GetNrows());
      gc2o= new TGraphErrors (cmp2x, cmp2obs, zero2, cmp2err);
    } else {
      gc2o= new TGraph (cmp2x, cmp2obs);
    }
    gc2o->SetLineWidth(2);
    gc2o->SetLineColor(kBlue);
    gc2o->SetMarkerColor(kBlue);
    gc2o->SetMarkerSize(0.8);
    gc2o->SetMarkerStyle(20);

    gc2c= new TGraph (cmp2x, cmp2c);
    gc2c->SetLineStyle(2);
    gc2c->SetLineWidth(2);
    gc2c->SetLineColor(6);
  }

  if (!c1) c1= new TCanvas("limitPlot","limitPlot",1050,750);
  gPad->SetLogx((logmu&4)?1:0);
  gPad->SetLogy((logmu&1)?1:0);
  gPad->SetTicks(0,1);
  gStyle->SetOptStat(0);

  TH2F* clplot= new TH2F ("limits", ";m_{H} [GeV];95% CL limit on #sigma/#sigma_{SM}", 50, mlo, mhi, 50, mulo, muhi);

  clplot->Draw();

  gb2->Draw("e3");
  gb1->Draw("e3");
  gc->Draw("l");
  if (go) go->Draw("l");
  gp->Draw("lp");
  if (compare) {
    gco->Draw("p");
    if (gcp1s) gcp1s->Draw("l");
    if (gcc) gcc->Draw("l");
    if (gcm1s) gcm1s->Draw("l");
  }
  if (compare2) {
    gc2o->Draw("p");
    //    gc2c->Draw("l");
  }

  clplot->Draw("axis,same");

  TLegend* leg = new TLegend (0.68, 0.74, 0.93, 0.93);
  if (useCLs) {
    leg->AddEntry(gp,Form("Obs. Asym. CL_{s} %s",nam.Data()),"lp");
    leg->AddEntry(gc,Form("Exp. Asym. CL_{s} %s",nam.Data()),"l");
  } else {
    leg->AddEntry(gp,Form("Obs. PCL %s",nam.Data()),"lp");
    if (go) leg->AddEntry(go,Form("Obs. %s",nam.Data()),"l");
    leg->AddEntry(gc,Form("Exp. %s",nam.Data()),"l");
  }
  leg->AddEntry(gb1,"#pm 1#sigma","f");
  if (useCLs)
    leg->AddEntry(gb2,"#pm 2#sigma","f");
  else
    leg->AddEntry(gb2,"+ 2#sigma","f");
  if (compare) {
    TString compname;
    if (cmpname)
      compname = cmpname;
    else {
      compname = gSystem->BaseName (compare);
      if (compname.EndsWith(".root")) compname.Remove(compname.Length()-5,5);
    }
    leg->AddEntry(gco,Form("Obs. %s",compname.Data()),"p");
    if (haveCmpBands) {
      if (haveCmpBands==1)
        leg->AddEntry(gcc,Form("Exp. %s",compname.Data()),"l");
      else if (useCLs)
        leg->AddEntry(gcc,Form("Exp. %s (median,#pm 1#sigma)",compname.Data()),"l");
      else
        leg->AddEntry(gcc,Form("Exp. %s (median,#pm 1#sigma)",compname.Data()),"l");
    }
  }
  if (compare2) {
    TString compname;
    if (cmpname)
      compname = cmp2name;
    else {
      compname = gSystem->BaseName (compare2);
      if (compname.EndsWith(".root")) compname.Remove(compname.Length()-5,5);
    }
    leg->AddEntry(gc2o,Form("Obs. %s",compname.Data()),"p");
    //    leg->AddEntry(gc2c,Form("Exp. %s",compname.Data()),"l");
  }
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->Draw("same");

  if ((logmu&1) || (mulo<1 && muhi<100)) {
    TLine* line = new TLine (mlo, 1.0, mhi, 1.0);
    line->SetLineColor(4);
    line->SetLineStyle(3);
    line->Draw();
  }

  ATLASLabel(0.15,0.04,true,1);

  TString psnam;
  psnam.Form("H_%s%s%s.ps",(nam.Length()==0 ? "limit" : nam.Data()),(!haveCLs?"":useCLs?"_cls":"_pcl"),((logmu&1)?"":"_lin"));
  c1->Print (psnam);
}
