import ROOT, os, sys, re, inspect
from os.path import dirname, normpath

class ldcms_dummy: pass

def ld (opt="k"):
  #  Load classes required by ATLAS Higgs workspaces
  #  The possible options are:
  #     k : keep the shared library after the session end. (default)
  #     f : force recompilation.
  #     g : compile with debug symbol
  #     O : optimized the code (ignore if 'g' is specified)
  #     c : compile only, do not attempt to load the library.
  ADD_WSCLASSES= os.environ.get("ADD_WSCLASSES","")
  if ADD_WSCLASSES=="1" or ADD_WSCLASSES=="2":
    ROOT.gSystem.CompileMacro("RooBSplineBases.cxx",            opt)
    ROOT.gSystem.CompileMacro("RooBSpline.cxx",                 opt)
    ROOT.gSystem.CompileMacro("RooParamKeysPdf.cxx",            opt)
    ROOT.gSystem.CompileMacro("RooStarMomentMorph.cxx",         opt)
  if ADD_WSCLASSES=="2":
    ROOT.gSystem.CompileMacro("RooLnuqqHighmass.cxx",           opt)
    ROOT.gSystem.CompileMacro("RooLnuqqSignal.cxx",             opt)
    ROOT.gSystem.CompileMacro("RooScaleLOSM.C",                 opt)
  if ADD_WSCLASSES=="3":
    ROOT.gSystem.CompileMacro("ABWxG.cxx",                      opt)
  if ADD_WSCLASSES=="4":
    ROOT.gSystem.CompileMacro("Background.cxx",                 opt)

def ldcms (opt="k"):
  #  Load classes required by CMS Higgs workspaces
  if os.environ.get("CMSCODE","")=="": return

  cmsdir= normpath (dirname (inspect.getsourcefile(ldcms_dummy))) + "/cmscode/src/";
  ROOT.gSystem.AddIncludePath("-I"+normpath(cmsdir+"../interface"));
  ROOT.gSystem.SetBuildDir(".",1);  # keep shared libs in this dir

  ROOT.gSystem.Load("libSmatrix");
  ROOT.gSystem.CompileMacro (cmsdir+"HGGRooPdfs.cc",            opt);
  ROOT.gSystem.CompileMacro (cmsdir+"HZZ2L2QRooPdfs.cc",        opt);
  ROOT.gSystem.CompileMacro (cmsdir+"HZZ4LRooPdfs.cc",          opt);
  ROOT.gSystem.CompileMacro (cmsdir+"ProcessNormalization.cc",  opt);
  ROOT.gSystem.CompileMacro (cmsdir+"RooBernsteinFast.cc",      opt);
  ROOT.gSystem.CompileMacro (cmsdir+"RooMultiPdf.cxx",          opt);
  ROOT.gSystem.CompileMacro (cmsdir+"RooSpline1D.cc",           opt);
  ROOT.gSystem.CompileMacro (cmsdir+"AsymPow.cc",               opt);

  # Tell ACLiC to use C++11, required by the following classes
  makeso= ROOT.gSystem.GetMakeSharedLib();
  ROOT.gSystem.SetMakeSharedLib(re.sub(r" -c ", " -std=c++0x -c ", makeso));

  ROOT.gSystem.CompileMacro (cmsdir+"SimpleCacheSentry.cc",     opt);
  ROOT.gSystem.CompileMacro (cmsdir+"FastTemplate.cc",          opt);
  ROOT.gSystem.CompileMacro (cmsdir+"VerticalInterpHistPdf.cc", opt);

  ROOT.gSystem.SetMakeSharedLib(makeso);
  ROOT.gSystem.SetBuildDir("");

ld()
ldcms()
