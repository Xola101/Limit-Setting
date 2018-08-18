//
// ATLAS Style, based on a style file from BaBar
//

#ifndef __ATLASSTYLE_C
#define __ATLASSTYLE_C

#include "AtlasStyle.h"

#include "TROOT.h"
#include "TStyle.h"

void SetAtlasStyle (bool force)
{
  static TStyle* atlasStyle = 0;
  // std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  if ( atlasStyle==0 ) atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  if (force) gROOT->ForceStyle();
}

TStyle* AtlasStyle(const char* baseStyle)
{
  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");
  if (baseStyle) {
    TStyle* base = gROOT->GetStyle(baseStyle);
    if (base) {
      base->Copy(*atlasStyle);
      // Following missed by copy
      atlasStyle->SetLegendFillColor(0);
      atlasStyle->SetLegendFont(42);
    }
  }

  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);      // Modern default
  atlasStyle->SetFrameFillColor(icol);       // Modern default
  atlasStyle->SetCanvasBorderMode(icol);     // Modern default
  atlasStyle->SetCanvasColor(icol);          // Modern default
  atlasStyle->SetPadBorderMode(icol);        // Modern default
  atlasStyle->SetPadColor(icol);             // Modern default
  atlasStyle->SetStatColor(icol);            // Modern default
  //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);           // default

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.4);
  atlasStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);            // Modern default
  atlasStyle->SetLabelFont(font,"x");        // Modern default
  atlasStyle->SetTitleFont(font,"x");        // Modern default
  atlasStyle->SetLabelFont(font,"y");        // Modern default
  atlasStyle->SetTitleFont(font,"y");        // Modern default
  atlasStyle->SetLabelFont(font,"z");        // Modern default
  atlasStyle->SetTitleFont(font,"z");        // Modern default
  
  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2);
  atlasStyle->SetHistLineColor(1);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes  - default

  // get rid of X error bars 
  //atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  return atlasStyle;

}

#endif // __ATLASSTYLE_C
