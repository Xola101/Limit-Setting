#!/usr/bin/env python
# $Id: npnote.py 765371 2016-07-29 15:43:04Z adye $

__author__  = "Tim Adye"
__version__ = "$Revision: 765371 $"

import sys, getopt, re, math, tempfile, glob, optparse, array
from sys import argv, exit, stdout
from os.path import dirname, basename
from math import sqrt, pi

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions= True   # prevent ROOT screwing with out help text
# Load RooFit without any startup message and make RooAbsCollection iterable.
from PyRootUtils import Struct, TDirectoryGlobIter, runcmd, CompileCpp, copysign
from RooPyUtils  import Corr
from ROOT import gStyle, gPad, gSystem, TVectorD
from ROOT import RooFit, RooRealVar

# Load custom workspace classes
import ld

prog= basename(argv[0])
opt=None
interactive=False
verbose= 0
tmpdir= None
__version__ = re.sub(r"^.* ([\d.]+).*$",r"\1",__version__)


def parseArgs():
  global opt, verbose, interactive
  results_name_def1= "minos_results"
  results_name_def2= "fitres/data_stage1_try? fit_result"   # use fit_result if defined, else info from fitres dir

  parser= optparse.OptionParser(usage="%prog [OPTIONS] IN.root [IN.root...]", version="%prog "+__version__)
  parser.add_option ("-v", "--verbose",      help="verbose running",                              action="count")
  parser.add_option ("-i", "--interactive",  help="ROOT interactive mode",                        action="store_true")
  parser.add_option ("-o", "--outname",      help="Outfile file NAME -> NAME.pdf etc [default=%default]", metavar="NAME", default="np")
  parser.add_option ("-r", "--results_file", help="MINOS results files (space-separated wildcards)", metavar="FILE")
  parser.add_option ("-2", "--results_file2",help="comparison results files", metavar="FILE")
  parser.add_option ("-3", "--results_file3",help="overriding correlation matrix results files (rank by separate fits)", metavar="FILE")
  parser.add_option ("-t", "--txtfile",      help="text file with Stefan Gadatsch's results",     metavar="FILE")
  parser.add_option ("-C", "--poi",          help="POI name(s) and LaTeX [default=%default]",     metavar="NAME")
  parser.add_option ("-N", "--nxy",          help="Plots per page in OUT.pdf [default=%default]", metavar="NxM", default="2x2")
  parser.add_option ("-p", "--pagexy",       help="Plots per page in OUT_note.pdf [default=%default]", metavar="NxM", default="4x5")
  parser.add_option ("-R", "--results_name", help="Object name in results file [default='%s' or '%s' if --use_migrad_errs=2]" % (results_name_def1, results_name_def2), metavar="NAME")
  parser.add_option ("-e", "--maxEnt",       help="Number of ntuple entries to process", metavar="N", default=-1, type="int")
  parser.add_option ("-m", "--maxPlots",     help="Maximum number of plots to produce",  metavar="N", default=-1, type="int")
  parser.add_option (      "--plotMuHatPositive",  help="Plot NP for muhat>=0",                action="store_true")
  parser.add_option (      "--include_corr",       help="0=don't include correlations, 1=print corr, 2=sort by corr",
                                                   metavar="N", type="int", default=1)
  parser.add_option (      "--pullsPerPlot",       metavar="N", type="int", default= 20)
  parser.add_option (      "--use_migrad_errs",    help="0=MINOS only, 1=Migrad errors, 2=Migrad covarance matrix",
                                                   metavar="N", type="int", default=0)
  parser.add_option (      "--all_cond",           action="store_true")
  parser.add_option (      "--no_nt_data",         action="store_true")
  parser.add_option ("-c", "--cut", default="")
  parser.add_option (      "--sf",                 metavar="N", type="int", default=3)
  parser.add_option (      "--no_latex",           action="store_true")
  parser.add_option (      "--no_sign",            help="don't indicate anticorrelations in blue", action="store_true", default=False)
  parser.add_option (      "--sign_up",            help="impact plot sign shows NP +/-1sigma, rather than +/-Delta_mu", action="store_true", default=False)
  parser.add_option (      "--oldsign",            help="Old (bad) sign convention: +UP/HI -DO/LO", action="store_true", default=False)
  parser.add_option (      "--pubnames",           action="store_true", default=False)
  parser.add_option (      "--sumrest",            metavar="N", type="int", default=-1)
  parser.add_option (      "--setval",             help="INDEX=NAME+HI-LO")
  parser.add_option (      "--maxsys",             metavar="N", type="int", default=-1)
  parser.add_option (      "--atlasLabel",         action="store_true", default=False)
  parser.add_option (      "--preliminary",        action="store_true", default=False)
  parser.add_option (      "--internal",           action="store_true", default=False)
  parser.add_option ("-7", "--7TeV",               action="store_true", default=False, dest="e7TeV")
  parser.add_option ("-8", "--8TeV",               action="store_true", default=False, dest="e8TeV")

  opt, args= parser.parse_args()
  verbose= opt.verbose
  interactive= opt.interactive
  if not opt.results_name:
    if opt.use_migrad_errs==2: opt.results_name=results_name_def2
    else:                      opt.results_name=results_name_def1
  opt.atlasLabel= opt.atlasLabel or opt.preliminary or opt.internal
  if not opt.e7TeV and not opt.e8TeV:
    opt.e7TeV= True
    opt.e8TeV= True
  if len(args)==0 and not opt.results_file and not opt.txtfile:
    parser.print_help()
    exit(1)
  return args


def process (args,
             texfile_main=         "minos_table.tex",
             texfile_include=      "minos_include_%s.tex",
             texfile_include_corr= "minos_include_corr.tex",
             texfile_figs=         "minos_figs.tex",
             texfile_include_figs= "minos_include_figs.tex",
             texfile_include_pulls="minos_include_pulls.tex"):

  ROOT.TH1.AddDirectory(0)
  if not interactive: ROOT.gROOT.SetBatch(True)

  pll= None
  for fname in args:
    f= ROOT.TFile.Open (fname)
    if not f: return 1
    ROOT.gROOT.cd()
    pll= MergeTree (f, pll)
    if not pll: return 2
  if pll: print "Ntuple '%s' has %d entries" % (pll.GetName(), pll.GetEntries())

  poiName= []
  poiLatex= []
  if opt.poi:
    for s in opt.poi.split(","):
      m= re.search(r"^([^=]+)=(.*)$",s)
      if m:
        n,l= m.groups()
      else:
        n=l=s
        l= latexName(l)
      poiName.append(n)
      poiLatex.append(l)

  use_only_nt_data= not opt.results_file
  opt.haveimp= 0
  opt.havecont= 0
  if use_only_nt_data:
    opt.no_nt_data= False
    err, vars, corfit= fakeDataFit(pll)
  else:
    err, vars, corfit= readDataFit(poiName,poiLatex)
  if err: return err

  if len(poiName)==0:
    poiName=["mu"]
    poiLatex=["\\mu"]

  if not vars:
    print "No",opt.results_name
    return 4

  nvars= len(vars)
  opt.have_corfit= bool(corfit)
  if not corfit: corfit= [nvars*[0.0] for i in xrange(nvars)]

  if not opt.no_nt_data: addNtupleData (pll, vars, use_only_nt_data)

  ilist, grp= groupVars (vars, poiName)
  nvars= len(ilist)
  ncons= grp[0].next
  ipoi=  grp[3].first
  vars= [vars[i] for i in ilist]
  ncorfit= len(corfit)
  corfit= [[(i<ncorfit and j<ncorfit) and corfit[i][j] or 0.0 for j in ilist] for i in ilist]

  if opt.txtfile: addTextData(vars[:ncons], vars[ipoi])

  ctexfile= opt.outname+"_corr_table.tex"
  cor= getCorr (pll, vars, corfit, ipoi, ctexfile)
  if opt.have_corr: corsort= cor
  else:             corsort= corfit
  if not opt.have_corr and not opt.have_corfit: opt.include_corr= 0

  if   opt.haveimp:
    isort= sortVarsImp (vars, grp)
  elif opt.include_corr:
    isort= sortVars (vars, corsort, grp)
    if opt.include_corr>=2:
      vars= [vars[i] for i in isort]
      cor= [[cor[i][j] for j in isort] for i in isort]
      corfit= [[corfit[i][j] for j in isort] for i in isort]
      isort= [i for i in xrange(nvars)]  # no longer needs sorting
  else:
    isort= [i for i in xrange(nvars)]

  nplot= PlotDists (pll, vars, grp, poiName, cor[ipoi:], corfit[ipoi:])

  texfile_template= opt.outname + "_table_%s.tex"
  texname_figs=     opt.outname + "_figs.tex"
  if not opt.no_latex:
    for g in grp[:-1]:
      writeLatexTable (texfile_template % g.name, vars[g.first:g.next], poiLatex)

    ftex_template= open(texfile_figs).read()
    if not ftex_template: return 5
    ftex= open (texname_figs, "w")
    if opt.ntoys:
      npnames= [v.name for v in vars]
      for g in grp[:-1]:
        writeLatexFigs (ftex, [vars[i].name for i in xrange(g.first,g.next)], ftex_template, g.name, g.title)
    del ftex

  texname_pulls= opt.outname + "_pulls.tex"
  if ncons>0: PullPlot (vars[:ncons], isort, texname_pulls, poiLatex)

  if not opt.no_latex:
    err= RunLaTeX (texfile_main, texfile_include, texfile_include_corr, texfile_include_figs, texfile_include_pulls,
                   texfile_template, ctexfile, texname_figs, texname_pulls, tmpdir, "figures/"+opt.outname, opt.outname+"_note.pdf")
    if err: return err

    # tar eps files
    tarfile= opt.outname+".tar.gz"
    print "Create tar file",tarfile,"from",tmpdir+"/*.eps"
    err= runcmd ("""\
cd '%s'
tar zcf '%s/%s' *.eps
rm -f *.eps
rmdir '%s'
""" % (tmpdir, gSystem.WorkingDirectory(), tarfile, tmpdir))
    if err: return err

  return 0


def niceName(name):
  varsubs= [
     [r"^alpha_",     ""],
     [r"^ATLAS_Hgg_", ""],
     [r"^ATLAS_",     ""],
     [r"_CMS_COMB$",  ""],
     [r"_ATLAS_COMB$",""],
#     [r"^ATLAS_(.*)",                   r"\1 A"],
#     [r"^CMS_(.*)",                     r"\1 C"],
#     [r"^H(gg|ZZ|WW|tt|bb)_(.* [AC])$", r"\2H\1"],
#     [r"^ttH_H(bb|tt|gg)_(.* [AC])$",   r"\2ttH\1"],
#     [r"^ttH_(lep|l)_(.* [AC])$",       r"\2ttH\1"],
#     [r"^ttH_(.* [AC])$",               r"\1ttH"],
#     [r" (\S+)$",                       r" #bf{\1}"],
  ]
  if opt.pubnames:
    varsubs += [
      [r"^EM_ES_Z$",                     r"Z#rightarrowee calibration"],
      [r"^EM_Gain$",                     r"LAr electronics non-linearity"],
      [r"^EM_L1Gain$",                   r"#splitline{LAr cell non-linearity}{(layer 1)}"],
      [r"^EM_L2Gain$",                   r"#splitline{LAr cell non-linearity}{(layer 2)}"],
      [r"^EM_LatLeak$",                  r"Lateral shower shape"],
      [r"^EM_LatLeakConv$",              r"Lateral shower shape (conv)"],
      [r"^EM_LatLeakUnconv$",            r"Lateral shower shape (unconv)"],
      [r"^EM_LArCalib_Barrel$",          r"LAr layer calibration (barrel)"],
      [r"^EM_LArUnconvCalib_Barrel$",    r"#splitline{LAr syst on material}{after presampler (barrel)}"],
      [r"^EM_LArElecUnconv_Barrel$",     r"#splitline{LAr syst on material}{before presampler (barrel)}"],
      [r"^EM_MatID_1$",                  r"ID material model (|#eta| < 1.1)"],
      [r"^EM_mRes_CT$",                  r"#splitline{Z measurement of}{energy resolution}"],
      [r"^EM_mRes_MAT$",                 r"#splitline{Energy resolution}{from material}"],
      [r"^EM_mRes_ST$",                  r"#splitline{Calorimeter intrinsic}{energy resolution}"],
      [r"^EM_BCKG_MoriondCat1$",         r"#splitline{H#rightarrow#gamma#gamma background model}{(unconv central low p_{Tt})}"],
      [r"^EM_BCKG_MoriondCat3$",         r"#splitline{H#rightarrow#gamma#gamma background model}{(unconv rest low p_{Tt})}"],
      [r"^EM_PS_Barrel$",                r"#splitline{Presampler energy scale}{(barrel)}"],
      [r"^isEM$",                        r"Photon identification"],
      [r"^MU_MS$",                       r"Muon momentum scale"],
      [r"^QCDscale_ggH$",                r"QCD scale, gg#rightarrowH"],
      [r"^pdf_Higgs_ggH$",               r"PDF+#alpha_{S}, gg#rightarrowH"],
      [r"^norm_SF_llll_Z_2mu2e_(201.)$", r"#splitline{Reducible background norm,}{H#rightarrowZZ*#rightarrow2#mu2e, \1}"],
      [r"^norm_SF_llll_Z_4e_(201.)$",    r"#splitline{Reducible background norm,}{H#rightarrowZZ*#rightarrow4e, \1}"],
      [r"^norm_SF_llll_Zbb_4mu_(201.)$", r"#splitline{Reducible background norm,}{H#rightarrowZZ*#rightarrow4#mu, \1}"],
      [r"^norm_SF_H4l_Z_llee_(201.)$",    r"#splitline{Reducible #font[12]{ll}ee background}{normalisation, \1}"],
      [r"^norm_SF_H4l_Zbb_llmumu_(201.)$",r"#splitline{Reducible #font[12]{ll}#mu#mu background}{normalisation, \1}"],
      [r"^MU_EFF$",                      r"Muon reco/id efficiencies"],
      [r"^EL_(201.)_IDST_high$",         r"#splitline{Electron reco/id efficiency,}{E_{T}>20 GeV, \1}"],
      [r"^LUMI_(201.)$",                 r"Luminosity, \1"],
      [r"^EL_(201.)_ST_15$",             r"#splitline{Electron reco/id efficiency,}{uncorr, 15<E_{T}<20 GeV, \1}"],
      [r"^BR_VV$",                       r"H#rightarrowZZ* branching ratio"],
      [r"^EL_(201.)_REC_low$",           r"#splitline{Electron reco efficiency,}{corr, 7< E_{T}<15 GeV, \1}"],
      [r"^pdf_qq$",                      r"PDF+#alpha_{S}, qq#rightarrowZZ*"],
      [r"^QCDscale_VV$",                 r"QCD scale, qq#rightarrowZZ*"],
      [r"^H4l_EL_EFF_ISOIP$",            r"#splitline{Electron isolation / impact}{parameter cut efficiency}"],
    ]
  for s in varsubs:
    name= re.sub (s[0], s[1], name)
  return name


def varsubs(name,allvars):
  if name in allvars: return name
  varsubs= [
    [r"^ATLAS_Hbb_(.*)",                r"alpha_\1_bb_1112"            ],
    [r"^ATLAS_Hbb_SysJetEResol_",       r"ATLAS_Hbb_JER_"              ],
    [r"^ATLAS_Hbb_Sys",                 r"ATLAS_Hbb_"                  ],
    [r"^ATLAS_SF_H4l_([^_]+)_([^_]+)_", r"ATLAS_SF_H4l_\2_\1_"         ],
    [r"^QCDscale_Bkg_VV_Htt$",          r"QCDscale_Bkg_VV"             ],
    [r"^alpha_Systtbar3JNorm_bb1112$",  r"alpha_Systtbar3JNorm_bb_1112"],
    [r"^pdf_Bkg_bb_VH$",                r"pdf_Higgs_VH"                ],

    [r"^alpha_(.*)_bb_1112$",           r"ATLAS_Hbb_\1"                ],
    [r"^ATLAS_Hbb_JER_",                r"ATLAS_Hbb_SysJetEResol_"     ],
    [r"^ATLAS_Hbb_",                    r"ATLAS_Hbb_Sys"               ],
    [r"^QCDscale_Bkg_VV$",              r"QCDscale_Bkg_VV_Htt"         ],
    [r"^alpha_Systtbar3JNorm_bb_1112$", r"alpha_Systtbar3JNorm_bb1112" ],
    [r"^pdf_Higgs_VH$",                 r"pdf_Bkg_bb_VH"               ],
  ]
  n= name
  for s in varsubs:
    n, rep= re.subn (s[0], s[1], name)
    if rep and n in allvars: break
  return n

def latexName(name):
  latexNames= [
    [r"mu",          r"\\mu"],
    [r"mH",          r"m_{H}"],
    [r"mu_(.*)",     r"\\mu_{\1}"],
    [r"lambda_gamZ", r"\\lambda_{\\gamma Z}"],
    [r"lambda_tauZ", r"\\lambda_{\\tau Z}"],
    [r"kappa_(.*)",  r"\\kappa_{\1}"],
    [r"lambda_(.*)", r"\\lambda_{\1}"],
  ]
  for s in latexNames:
    sub, nsub= re.subn ("^"+s[0]+"$", s[1], name)
    if nsub: return sub
  return name

def readDataFit(poiName,poiLatex):
  if opt.poi: poiname=poiName[0]
  else:       poiname="mu"
  vars= []
  allvars={}
  corfit=None
  poiimplo= poiimphi= poicontlo= poiconthi= None
  poiimpbad= poicontbad= False
  CompileCpp("GlobIter.h")  # define ROOT.GlobIter.kIgnore in advance
  for isres2,resfile in enumerate([opt.results_file, opt.results_file2, opt.results_file3]):
    if resfile is None: continue
    if isres2==2: corfit=None  # if -3 specified, only use this one for corfit
    found=[]
    for spec in resfile.split(" "):
      found += glob.glob(spec)
    if not found: found= [spec]   # Let TFile.Open complain
    for fname in found:
      f= ROOT.TFile.Open (fname)
      if not f: return 1, None, None
      if not opt.poi:
        fitted_poi= f.Get("fitted_poi")
        if fitted_poi:
          pname=[]
          for var in fitted_poi:
            if not isinstance(var,ROOT.RooRealVar): continue
            pname.append(var.GetName())
          if len(pname)>0:
            if len(poiName)>0:
              if pname!=poiName:
                print "WARNING:",fitted_poi.ClassName()+"::"+fitted_poi.GetName(),"in file",fname,"has different POIs:",",".join(pname)
            else:
              print fitted_poi.ClassName()+"::"+fitted_poi.GetName(),"in file",fname,"has POIs:",",".join(pname)
              poiName.extend(pname)
              poiname=poiName[0]
              poiLatex.extend(pname)
              for p in range(len(poiLatex)):
                poiLatex[p]= latexName(poiLatex[p])
      rname=None
      for rname in (name for spec in opt.results_name.split(" ") for name in TDirectoryGlobIter(f,spec,ROOT.GlobIter.kIgnore)): pass  # pick last one
      if rname: res= f.Get(rname)
      else:     res= None
      if not res:
        print "Could not find",opt.results_name,"in file",fname
      else:
        print res.ClassName()+"::"+res.GetName(),"in file",fname
        if isinstance(res,ROOT.RooFitResult):
          if corfit:
            print "Multiple correlation matrix results"
          else:
            M= res.correlationMatrix()
            corfit= [[M(i,j) for j in xrange(M.GetNcols())] for i in xrange(M.GetNrows())]
            del M
          res= ROOT.RooArgSet(res.floatParsFinal())
        if isres2==2: continue
        for var in res:
          if not isinstance(var,ROOT.RooRealVar): continue
          v=varInfo(var,isres2)
          name= v.name
          if isres2: name= varsubs(name,allvars)
          if name in allvars:
            if isres2 and name != v.name: print "Rename",v.name,"to",name
            a= vars[allvars[name]]
            if isres2:
              if v.errqual2>a.errqual2:
                for o in "var2 val2 err2 errlo2 errhi2 errqual2 errtype2 symm2 nf2 haveval2".split(): a[o]= v[o]
              elif v.errqual2==a.errqual2 and v.errqual2>0:
                print "Multiple MINOS results for",
                a.var2.Print()
                print "                       and",
                v.var2.Print()
            else:
              if v.errqual>a.errqual:
                vars[allvars[v.name]]= v
              elif v.errqual==a.errqual and v.errqual>0:
                print "Multiple MINOS results for",
                a.var.Print()
                print "                       and",
                v.var.Print()
          else:
            allvars[v.name]=len(vars)
            vars.append(v)

      impname= "impact_"+poiname
      imp= f.Get(impname)
      if not imp:
        print "No '%s' object in file %s" % (impname, fname)
      else:
        print imp.ClassName()+"::"+imp.GetName(),"in file",fname
        for var in imp:
          if not isinstance(var,ROOT.RooRealVar): continue
          name= var.GetName()
          if isres2: name= varsubs(name,allvars)
          if name == poiname:
            poiimplo= var.getAsymErrorLo()
            poiimphi= var.getAsymErrorHi()
            if poiimplo>=0 or poiimphi<=0:
              print "bad POI %s errors: %+g %+g - use %g" % (poiname, poiimphi, poiimplo, var.getError())
              poiimphi= var.getError()
              poiimplo= -poiimphi
              poiimpbad= True
            continue
          if name not in allvars:
            print "No variable",name
            continue
          a= vars[allvars[name]]
          if a.haveimp: continue
          anti= var.getAttribute("anticorrelated")
          ss=   var.getAttribute("same-sign")
          a.impflag= ""
          if anti: a.impflag += "A"
          if ss:   a.impflag += "S"
          if a.impflag=="": a.impflag="-"   # so we can keep still use fields in printout
          a.imphi= abs(var.getAsymErrorHi())
          a.implo= abs(var.getAsymErrorLo())
          if anti: (a.implo, a.imphi)= (a.imphi, -a.implo)
          if ss:    a.implo=  copysign (a.implo,  a.imphi)
          else:     a.implo=  copysign (a.implo, -a.imphi)
          a.impvar= var
          a.haveimp= 1
          opt.haveimp= 1

      contname= "contrib_"+poiname
      cont= f.Get(contname)
      if not cont:
        print "No '%s' object in file %s" % (contname, fname)
      else:
        print cont.ClassName()+"::"+cont.GetName(),"in file",fname
        for var in cont:
          if not isinstance(var,ROOT.RooRealVar): continue
          name= var.GetName()
          if isres2: name= varsubs(name,allvars)
          if name == poiname:
            poicontlo= var.getAsymErrorLo()
            poiconthi= var.getAsymErrorHi()
            if poicontlo>=0 or poiconthi<=0:
              print "bad POI %s errors: %+g %+g - use %g" % (poiname, poiconthi, poicontlo, var.getError())
              poiconthi= var.getError()
              poicontlo= -poiimphi
              poicontbad= True
            continue
          if name not in allvars:
            print "No variable",name
            continue
          a= vars[allvars[name]]
          if a.havecont: continue
          a.contlo= var.getAsymErrorLo()
          a.conthi= var.getAsymErrorHi()
          a.havecont= 1
          opt.havecont= 1

  if opt.haveimp:
#    topvars= """ATLAS_EM_Gain ATLAS_EM_LatLeak ATLAS_EM_LArCalib_Barrel ATLAS_EM_ES_Z ATLAS_EM_LArElecUnconv_Barrel
#                ATLAS_EM_MatID_1 ATLAS_EM_mRes_CT ATLAS_EM_mRes_MAT ATLAS_EM_BCKG_MoriondCat3 ATLAS_EM_mRes_ST
#                QCDscale_ggH pdf_Higgs_ggH ATLAS_MU_MS ATLAS_EM_BCKG_MoriondCat1""".split()
    if poiimplo is None:
      if poiname in allvars:
        poi= vars[allvars[poiname]]
        poiimplo= -poi.errlo
        poiimphi=  poi.errhi
        print "no POI %s in impact_%s, using %g %+g %+g" % (poiname, poiname, poi.val, poiimphi, poiimplo)
      else:
        print "POI %s not found" % poiname
        poiimplo= poiimphi= None
        opt.haveimp= 0
    elif poiimpbad:
      print "POI %s Migrad errors: +/- %g" % (poiname, poiimphi)
      opt.haveimp=2
    else:
      print "POI %s errors: %+g %+g" % (poiname, poiimphi, poiimplo)
    totimphi= totimplo= totimp= totimpmx= 0.0
    print "%-45s %9s %9s %-2s %9s %9s %9s %9s %9s %9s" % tuple("NP VARHI VARLO AS UP DO HI LO AV MX".split())
    for v in vars:
      if not v.haveimp: continue
      poi= vars[allvars[poiname]]
      if v.imphi>v.implo: imphi, implo= v.imphi, v.implo
      else:               imphi, implo= v.implo, v.imphi
#      imphi *= 100.0/poi.val
#      implo *= 100.0/poi.val
      imp=        (abs(implo)+ abs(imphi))/2.0
      impmx = max (abs(implo), abs(imphi))
#      if v.name in topvars:
      totimphi += imphi*imphi
      totimplo += implo*implo
      totimp   += imp  *imp
      totimpmx += impmx*impmx
      print "%-45s %+9.6f %+9.6f %-2s %+9.6f %+9.6f %+9.6f %+9.6f %9.6f %9.6f" % (v.name, v.impvar.getAsymErrorHi(), v.impvar.getAsymErrorLo(), v.impflag, v.imphi, v.implo, imphi, implo, imp, impmx)
      if poiimplo is None:
        del v.implo, v.imphi
        v.haveimp= 0
      else:
        if opt.oldsign:
          v.implo /= abs (poiimplo)
          v.imphi /= abs (poiimphi)
        else:
          if v.implo>=0.0: fimplo= poiimphi
          else:            fimplo= poiimplo
          if v.imphi>=0.0: fimphi= poiimphi
          else:            fimphi= poiimplo
          if fimplo==0.0: v.implo= 0.0
          else:           v.implo /= abs(fimplo)
          if fimphi==0.0: v.imphi= 0.0
          else:           v.imphi /= abs(fimphi)
        v.impmx = max (abs(v.implo), abs(v.imphi))
        if poiimpbad: v.haveimp= 2
    print "%-45s %+9.6f %+9.6f %-5s %6s %9s %+9.6f %+9.6f %9.6f %9.6f" % ("Total", poiimphi, poiimplo, "(POI)", "-", "-", sqrt(totimphi), -sqrt(totimplo), sqrt(totimp), sqrt(totimpmx))
  if not opt.haveimp: opt.no_sign= True

  if opt.havecont:
    if poicontlo is None:
      if poiname in allvars:
        poi= vars[allvars[poiname]]
        poicontlo= -poi.errlo
        poiconthi=  poi.errhi
        print "no POI %s in contrib_%s, using %g %+g %+g" % (poiname, poiname, poi.val, poiconthi, poicontlo)
      else:
        print "POI %s not found" % poiname
        poicontlo= poiconthi= None
        opt.havecont= 0
    elif poicontbad:
      print "POI %s Migrad errors: +/- %g" % (poiname, poiconthi)
      opt.havecont= 2
    else:
      print "POI %s errors: %+g %+g" % (poiname, poiconthi, poicontlo)
    for v in vars:
      if not v.havecont: continue
      if poicontlo is None:
        del v.contlo, v.conthi
        v.havecont= 0
      else:
        print "%s = %g %+g %+g" % (v.name, v.val, v.conthi, v.contlo),
        if v.contlo > poicontlo: v.contlo= sqrt (1.0 - (v.contlo/poicontlo)**2)
        else:                    v.contlo= 0.0
        if v.conthi < poiconthi: v.conthi= sqrt (1.0 - (v.conthi/poiconthi)**2)
        else:                    v.conthi= 0.0
        v.contmx = max (v.contlo, v.conthi)
        print "-> %g %g" % (v.conthi, v.contlo)
        if poicontbad: v.havecont= 2

  varlist= ROOT.RooArgList()
  for v in vars:
    if v.haveval and (not opt.results_file2 or v.haveval2): varlist.add(v.var)
  varlist.sort()
  print "All variables:"
  varlist.Print("v")

  if opt.results_file2:
    varlist= ROOT.RooArgList()
    for v in vars:
      if v.haveval and not v.haveval2: varlist.add(v.var)
    if varlist.getSize() > 0:
      print "Only in",opt.results_file+":-"
      varlist.sort()
      varlist.Print("v")

    varlist= ROOT.RooArgList()
    for v in vars:
      if not v.haveval and v.haveval2: varlist.add(v.var2)
    if varlist.getSize() > 0:
      print "Only in",opt.results_file2+":-"
      varlist.sort()
      varlist.Print("v")

  return 0, vars, corfit

def fakeDataFit(pll):
  # Fill with dummy vars made up from ntuple
  skip_var= r"^(weight|is_data|is_alt|sample|have_uncond|have_cond|mu_init|.*_cond|.*_err|globObs_.*)$"
  vars= []
  for leaf in pll.GetListOfLeaves():
    name= leaf.GetName()
    if re.match(skip_var,name): continue
    if not pll.GetBranch(name+"_err"): continue
    vars.append (varInfo (ROOT.RooRealVar (name, name, -999.0)))
  return 0, vars, None


def addNtupleData (pll, vars, use_only_nt_data):
  if not pll: return
  nent= pll.GetEntries()
  if opt.maxEnt>=0 and nent>opt.maxEnt: nent= opt.maxEnt
  pll.SetBranchStatus ("*", 0)
  pll.SetBranchStatus ("is_data", 1)
  for idata in xrange(0,nent):
    pll.GetEntry(idata)
    if pll.is_data: break
  pll.SetBranchStatus ("*", 1)
  if not pll.is_data:
    print "No data entry in ntuple"
    return
  pll.GetEntry(idata)
  suff=""
  if pll.FindLeaf("have_cond"):
    suff= "_cond"  # use conditional fit if we ran it
  else:
    print "No conditional fit in ntuple for data - use unconditional fit"
  for v in vars:
    cname= v.name+suff
    leaf= pll.FindLeaf(cname)
    if not leaf: continue
    v.val= leaf.GetValue()
    v.var.setVal(v.val)
    if use_only_nt_data and opt.use_migrad_errs:
      eleaf= pll.FindLeaf(cname+"_err")
      if eleaf:
        v.err= eleaf.GetValue()
        v.errtype="*"
        v.errlo= -v.err
        v.errhi=  v.err
        v.errqual=0
        v.var.setError(v.err)
        v.var.setConstant(0)
    v.have_data= True


def addTextData (vars, poi):
  varmap={}
  txttotal= None
  for v in vars: varmap[v.name]= v
  f= open(opt.txtfile)
  for il,s in enumerate(f):
    a= s.rstrip().split(" ")
    if len(a)==5 and a[0]=="total":
      if txttotal: print "%s:%d: duplicate: %s" % (opt.txtfile,il,s)
      else:        txttotal= [float(f) for f in a[2:5]]
    elif len(a)==8 or len(a)==10:
      v= varmap.get(a[0])
      if not v:
        print "%s:%d: unknown variable: %s" % (opt.txtfile,il,s)
      elif v.havetxt:
        print "%s:%d: duplicate: %s" % (opt.txtfile,il,s)
      else:
        v.havetxt= True
        v.txtvec= [float(f) for f in a[2:8]]
    else:
      print "%s:%d: bad format: %s" % (opt.txtfile,il,s)

  for v in vars:
    if not v.havetxt:
      print "%s: missing variable '%s'" % (opt.txtfile, v.name)

  if not txttotal:
    txttotal= [poi.val, -poi.errlo, poi.errhi]
    print "%s: no total error given, using %g +%g -%g" % (opt.txtfile, txttotal[0], txttotal[1], txttotal[2])

  totval, toterrhi, toterrlo= txttotal
  for v in vars:
    if not v.havetxt: continue
    mu, muhi, mulo= v.txtvec[3:6]
    if muhi>=mu:
      v.fiterrmx= max (abs((mu-mulo)/toterrlo), abs((muhi-mu)/toterrhi))
      v.fiterrlo= abs((mu-mulo)/toterrlo)
      v.fiterrhi= abs((muhi-mu)/toterrhi)
    else:
      v.fiterrmx= max (abs((mu-muhi)/toterrlo), abs((mulo-mu)/toterrhi))
      v.fiterrlo= abs((mu-muhi)/toterrlo)
      v.fiterrhi= abs((mulo-mu)/toterrhi)
    del v.txtvec

def groupVars (vars, poiName):
  # Sort into groups
  unconstrained= r"^(ATLAS_norm_(?!SF_).*|ATLAS_norm_SF_ZH.*|ATLAS_LUMI_.*_lvlv.*|a[1-4]_.*|a[12]EP_.*|p[0-3]_.*|atlas_nbkg_.*|slope_.*|ATLAS_Hbb_norm_.*|ATLAS_SF_HWW_WW0j_WW_.*|ATLAS_SF_HWW_WW1j_.*|ATLAS_SF_HWW_WW2j_.*|ATLAS_PM_EFF_f_recoil_.*|PM_EFF_f_recoil_.*|ATLAS_Systtbar3JNorm_7TeV.*|CMS_.*2j_stat_CMS_COMB|[sS]cale_.*|unconstr?_.*|u_.*|CMS_.*_SHAPE.*|mu|mu_.*|muVH|XSRatio|SigXsecOverSM|mH|m_4l|dM|mHiggs|DeltaM|kappa_.*|lambda_.*|ATLAS_qqqq_bkg_p[23])$"
  gammaStat= r"^(gamma_stat_.*|ATLAS_.*_gamma_stat_.*|CMS_.*_stat|CMS_.*_stat_CMS_COMB)$"

  grp= [Struct(i=[], name="constrained",   title="constrained"),                # 0
        Struct(i=[], name="unconstrained", title="unconstrained"),              # 1
        Struct(i=[], name="gamma",         title="histogram bin statistics"),   # 2
        Struct(i=[], name="POI",           title="POI")]                        # 3
  for v in vars:
    if            poiName.count(v.name): v.grp= 3
    elif re.match(unconstrained,v.name): v.grp= 1
    elif re.match(gammaStat,    v.name): v.grp= 2
    else:                                v.grp= 0
  for i,v in enumerate(vars):
    grp[v.grp].i.append(i)
  if len(grp[3].i)==0:
    grp[3].i= range(len(vars),len(vars)+len(poiName))
    for n in poiName:
      v= varInfo(RooRealVar(n,n,-999.0))
      v.grp= 3
      vars.append(v)
  if opt.maxPlots>=0:
    m= opt.maxPlots-len(grp[3].i)
    for g in grp[:-1]:
      if m<len(g.i): del g.i[max(0,m):]
      m -= len(g.i)
  ilist= [i for g in grp for i in g.i]
  j=0
  for g in grp:
    g.first= j
    j += len(g.i)
    del g.i
    g.next=  j
  stdout.write ("%d variables: %s: " % (len(ilist),", ".join(["%d %s"%(g.next-g.first,g.name) for g in grp])))
  for i in ilist[grp[3].first:grp[3].next]: vars[i].var.Print()
  return ilist, grp

def getCorr (pll, vars, corfit, ipoi, ctexfile):
  # NP-NP correlations
  if not opt.no_latex:
    ctex= open (ctexfile, "w")  # even if !include_corr, still make file to keep LaTeX happy
  npnames= [v.name for v in vars]
  cor= None
  if opt.include_corr and pll:
    cor= Corr (pll, npnames, opt.cut, opt.maxEnt)
  opt.have_corr= bool(cor)
  if not cor:
    nvars= len(vars)
    cor= [nvars*[0.0] for i in xrange(nvars)]
  if opt.have_corr or opt.have_corfit:
    print "%-50s %-50s" % ("Nuisance Parameter 1", "Nuisance Parameter 2"),
    if opt.have_corr:   print "%10s" % "rho(P1,P2)",
    if opt.have_corfit: print "%10s" % "MIGRAD rho",
    print
    if not opt.no_latex:
      print >>ctex, "%-50s& %-50s" % ("% Nuisance Parameter 1", "% Nuisance Parameter 2"),
      if opt.have_corr:   print >>ctex, "& %16s" % "$\\rho(P_1,P_2) \\%$",
      if opt.have_corfit: print >>ctex, "& %16s" % "$\\rho(P_1,P_2) \\%$",
      print >>ctex, "\\\\ \\hline"
    for i in xrange(0,ipoi):  # skip mu
      for j in xrange(i+1,ipoi):
        rho=    cor[i][j]
        rhofit= corfit[i][j]
        if abs(rho)>0.2 or abs(rhofit)>0.2:
          name1= npnames[i]
          name2= npnames[j]
          print "%-50s %-50s" % (name1, name2),
          if opt.have_corr:   print "%10.2f" % (100.0*rho),
          if opt.have_corfit: print "%10.2f" % (100.0*rhofit),
          print
          if not opt.no_latex:
            print >>ctex, "%-50s& %-50s" % (re.sub(r"_","\\_",name1), re.sub(r"_","\\_",name2)),
            if opt.have_corr:   print >>ctex, "& %16.2f" % (100.0*rho),
            if opt.have_corfit: print >>ctex, "& %16.2f" % (100.0*rhofit),
            print >>ctex, "\\\\ \\hline"
    print
  return cor


def sortVars (vars, cor, grp, txttotal=None):
  ipoi= grp[3].first
  xsort= [(-abs(cor[ipoi][i]),niceName(v.name).lower()) for i,v in enumerate(vars)]
  if opt.txtfile:
    for i,v in enumerate(vars):
      if v.havetxt:
        if v.fiterrmx>-xsort[i][0]: xsort[i]= (-v.fiterrmx,xsort[i][1])
  ixsort= list(enumerate(xsort))
  return [m[0] for g in grp for m in sorted (ixsort[g.first:g.next], key=lambda m:m[1])]

def sortVarsImp (vars, grp):
  ipoi= grp[3].first
  ixsort= []
  for i,v in enumerate(vars):
    if v.haveimp: ixsort.append((i,(-abs(v.impmx),niceName(v.name))))
    else:         ixsort.append((i,(0.0,niceName(v.name))))
  return [m[0] for g in grp for m in sorted (ixsort[g.first:g.next], key=lambda m:m[1])]



def PlotDists (pll, vars, grp, poiName, poicor, poicorfit):
  nx, ny= [int(n) for n in re.match(r"^(\d+)x(\d+)$",opt.nxy).groups()]

  gStyle.SetOptStat(0)
  canvas= ROOT.TCanvas (prog, prog, 700, 700)
  canvas.SetFillStyle(0);
  plotFile= opt.outname+".pdf"
  PlotFile(canvas,plotFile)

  havetoys= pll and pll.GetEntries()>0
  print "%-50s %10s%s %10s %10s%1s %10s" % ("Nuisance Parameter", "Data fit", " ", "MINOS -err", "MINOS +err", " ", "MIGRAD err"),
  if havetoys: print "%10s %10s" % ("Toy mean", "Toy sigma"),
  if opt.have_corr:
    for p in poiName:
      print "%10s" % ("rho(P,%s)%%" % p),
  if opt.have_corfit:
    for c in poicorfit:
      print "%10s" % "MIGRAD rho",
  if opt.haveimp:
    print "%21s" % (poiName[0]+"-impact (-/+)"),
  if opt.havecont:
    print "%21s" % (poiName[0]+"-contrib (-/+)"),
  print

  nvars= len(vars)
  page= []
  saved= []
  nplot= 0
  maxtoys= 0
  lastgrp= None

  for ivar,v in enumerate(vars):

    if v.grp != lastgrp:
      print "==========",grp[v.grp].title,"=========="
      lastgrp= v.grp

    varname= v.name
    if opt.all_cond: varname += "_cond"
    ntoys, h, draw= PlotNP (pll, varname, v.name, v.val, v.errlo, v.errhi, opt.cut, opt.maxEnt, opt.all_cond)
    if ntoys>maxtoys: maxtoys= ntoys

    toymean=0.0
    toyrms=0.0
    if ntoys:
      saved.append(h)
      nplot += 1
      toymean= h.GetMean()
      toyrms=  h.GetRMS()

      if not opt.no_latex: PlotEps (v.name, draw)

      f= "%%10.%df"%v.nf
      nw= max (len((f%toymean).lstrip()), len((f%toyrms).lstrip()), len((f%v.val).lstrip()))
      if v.symm:
        nw= max (nw, len((f%v.err).lstrip()))
      else:
        f= "%%+10.%df"%v.nf
        nw= max (nw, len((f%v.errhi).lstrip()), len((f%v.errlo).lstrip()))

      draw.append (DrawLatex   (0.89, 0.87, ("Toy mean  %%%d.%df"          % (nw,v.nf)) % toymean))
      draw.append (DrawLatex   (0.89, 0.84, ("Toy #sigma         %%%d.%df" % (nw,v.nf)) % toyrms))
      draw.append (DrawLatex   (0.89, 0.81, ("Data fit      %%%d.%df"      % (nw,v.nf)) % v.val))
      if   v.errtype=="C": pass
      elif v.errtype=="*":
        draw.append (DrawLatex (0.89, 0.78, ("Migrad err %%%d.%df"         % (nw,v.nf)) % v.err))
      elif v.errtype=="+":
        draw.append (DrawLatex (0.89, 0.78, ("MINOS err %%+%d.%df"         % (nw,v.nf)) % v.err))
      elif v.errtype=="-":
        draw.append (DrawLatex (0.89, 0.78, ("MINOS err %%+%d.%df"         % (nw,v.nf)) % -v.err))
      elif v.symm:
        draw.append (DrawLatex (0.89, 0.78, ("MINOS err %%%d.%df"          % (nw,v.nf)) % v.err))
      else:
        draw.append (DrawLatex (0.89, 0.78, ("MINOS err %%+%d.%df"         % (nw,v.nf)) % v.errhi))
        draw.append (DrawLatex (0.89, 0.75, (          "%%+%d.%df"         % (nw,v.nf)) % v.errlo))

      page.append(draw)
      if len(page)==nx*ny:
        PlotPdf (canvas, plotFile, page, nx, ny)
        page= []

    if opt.no_nt_data or v.have_data: cc=" "
    else:                             cc="*"
    print ("%%-50s %%10.%df%%c %%10.%df %%10.%df%%c %%10.%df" % (v.nf,v.nf,v.nf,v.nf)) % (v.name, v.val, cc, v.errlo, v.errhi, v.errtype, v.var.getError()),
    if ntoys:      print ("%%10.%df %%10.%df" % (v.nf,v.nf)) % (toymean, toyrms),
    elif havetoys: print "%10s %10s" % ("",""),
    if opt.have_corr:
      for c in poicor:
        print "%10.2f" % (100.0*c[ivar]),
    if opt.have_corfit:
      for c in poicorfit:
        print "%10.2f" % (100.0*c[ivar]),
    if opt.haveimp:
      if v.haveimp:
        print "%+10.2f %+10.2f" % (100.0*v.implo, 100.0*v.imphi),
        if not opt.havecont: print v.impflag,
    if opt.havecont:
      if v.havecont:
        print "%+10.2f %+10.2f" % (100.0*v.contlo, -100.0*v.conthi),
    print

    v.update (toymean=toymean, toyrms=toyrms, ntoys=ntoys,
              poicor=   [c[ivar] for c in poicor],
              poicorfit=[c[ivar] for c in poicorfit])

  if len(page)>0: PlotPdf (canvas, plotFile, page, nx, ny)
  PlotFile(canvas,plotFile,True)
  opt.ntoys= maxtoys

  if saved:
    resultFileName= opt.outname+".root"
    print "Write %d histograms to '%s'" % (len(saved), resultFileName)
    fileOut= ROOT.TFile.Open (resultFileName, "recreate")
    for h in saved: fileOut.WriteTObject (h)
    del fileOut

  return nplot


def writeLatexTable (fname, vars, poiLatex):
  tex= open (fname, "w")
  print >>tex, "%-50s&%10s &%30s &%10s &%12s" % ("% Nuisance Parameter", "Data fit", "Error", "Toy mean", "Toy $\\sigma$"),
  for p in poiLatex:
    if opt.include_corr: print >>tex, " &%16s" % ("$\\rho(%s,P) \\%%$" % p),
    else:                print >>tex, " &",
# if opt.haveimp or opt.havecont:
#   print >> tex, " &%21s" % ("$\\Delta\\hat{%s} / \\Delta\\hat{%s}_{\\mathrm{tot}}" % (poiLatex[0], poiLatex[0])),
  print >>tex, " \\\\ \\hline"

  for v in vars:
    name= re.sub(r"_","\\_",v.name)
    print >>tex, ("%%-50s&%%10.%df &"%v.nf) % (name, v.val),
    if not opt.use_migrad_errs and (v.err<=0.0 or v.errtype=='*'):
      print >>tex, "%28s   &" % "",
    elif v.symm:
      print >>tex, ("%%28.%df   &" % v.nf) % v.err,
    else:
      print >>tex, "%30s &" % (("$_{-%%.%df}^{+%%.%df}$" % (v.nf,v.nf)) % (abs(v.errlo), v.errhi)),
    if v.ntoys:
      print >>tex, "%%10.%df &%%12.%df" % (v.nf,v.nf) % (v.toymean, v.toyrms),
    else:
      print >>tex, "%10s &%10s" % ("",""),
    if opt.have_corr: pc= v.poicor
    else:             pc= v.poicorfit
    for c in pc:
      if opt.include_corr: print >>tex, " &%16.2f" % (100.0*c),
      else:                print >>tex, " &",
    print >>tex, " \\\\ \\hline"


def writeLatexFigs (tex, names, ftex_template, varType, varTypeName):
  pagex,pagey= [int(n) for n in re.match(r"^(\d+)x(\d+)$",opt.pagexy).groups()]
  npage=pagex*pagey
  nnames= len(names)
  ipage=0
  init=True
  data=""
  for j,name in enumerate(names):
    i= j+1
    if init:
      data += "\\mbox{%\n"
      init=False
    data += "\\includegraphics[width=%.3f\\textwidth]{figures/%s/%s.eps}%%\n" % (1.0/pagex, opt.outname, name)
    if not (i%pagex==0 or i==nnames): continue
    init=True
    if not (i%npage==0 or i==nnames):
      data += "}\\\\[0.5ex]\n"
    else:
      data += "}%\n"
      lab= ""
      if ipage==0:  lab += "\n\\label{fig:np_%s}"     % varType
      if i==nnames: lab += "\n\\label{fig:np_%s_end}" % varType
      print >>tex, ftex_template % (data, varTypeName, lab),
      data= ""
      ipage += 1


def PullPlot (vars, isort, texname_pulls, poiLatex, txttotal=None, texfile_pulls= "minos_pulls.tex"):

  pullFile= opt.outname+"_pulls.pdf"
  canvas= ROOT.TCanvas (prog, prog, 700, 800)
  canvas.SetFillStyle(0);
  PlotFile(canvas,pullFile)

  if opt.haveimp==2: print "WARNING: no MINOS errors on %s, used in relative impact" % poiLatex[0]
  if   opt.haveimp:      showcorr=1
  elif opt.include_corr: showcorr=2
  else:                  showcorr=0
  maxname= max([len(niceName(v.name)) for v in vars])

  if opt.pullsPerPlot>=15: nlines=5
  else:                    nlines=1

  if   opt.pullsPerPlot>15: bmargin= tmargin= 1.0
  elif opt.pullsPerPlot>8:  bmargin= tmargin= 0.4
  else:                     bmargin= tmargin= 0.25
  if opt.atlasLabel: tmargin += 1.0

  canvas.SetLogy(0)
  if   maxname>=45: left=0.45
  elif maxname>=35: left=0.35
  else:             left=0.31
  canvas.SetLeftMargin(left)
  canvas.SetRightMargin(0.03)
  top= 0.06
  canvas.SetTopMargin (top)

  if not opt.no_latex:
    ptex_template= open(texfile_pulls).read()
    if not ptex_template: return 5
    ptex= open (texname_pulls, "w")

  # Convert LaTeX label into a TLatex label "#Delta#hat{POI}(#theta=#pm1#sigma) / #sigma_{tot}(POI)", tidying sub/superscripts
  poilab= re.sub(r"\\","#",poiLatex[0])
  poilabLines=1
  if   opt.haveimp:
    poilab1, sub1=   re.subn (r"^([^_^]+)([_^])([^{])$",r"#hat{\1}\2{\3}", poilab)
    if not sub1:
      poilab1, sub1= re.subn (r"^([^_^]+)([_^])(\{.*\})$",r"#hat{\1}\2\3", poilab)
    if not sub1:
      poilab1= "#hat{%s}" % poilab
    poilab= "#frac{#Delta%s(#theta=#pm1#sigma)}{#sigma_{tot}(%s)}" % (poilab1, poilab)
    poilabLines=2
  elif opt.include_corr:
    poilab= "#rho(%s,P)" % poilab
  poilab += " [%]"
  if not opt.no_sign: poilabLines += 1
  if showcorr and not opt.txtfile: print "Yellow label:", poilab

  pullmax=2.0
  for v in vars:
    if v.val!=-999.0:
      if v.val+v.errlo      < -pullmax: pullmax= -1.01*(v.val+v.errlo)
      if v.val+v.errhi      >  pullmax: pullmax=  1.01*(v.val+v.errhi)
    if v.ntoys:
      if v.toymean-v.toyrms < -pullmax: pullmax= -1.01*(v.toymean-v.toyrms)
      if v.toymean+v.toyrms >  pullmax: pullmax=  1.01*(v.toymean+v.toyrms)

  ncons= len(vars)
  if   opt.haveimp:     cormax= max (0.1, abs(vars[isort[0]].impmx))
  elif opt.have_corr:   cormax= max (0.1, abs(vars[isort[0]].poicor[0]))
  else:                 cormax= max (0.1, abs(vars[isort[0]].poicorfit[0]))
  if opt.txtfile:
    for v in vars:
      if v.havetxt and v.fiterrmx>cormax: cormax= v.fiterrmx

  if opt.sumrest>=0:
    resthi= restlo= 0.0
    for i in range(opt.sumrest,ncons):
      v=vars[isort[i]]
      if v.imphi>=0.0: resthi += v.imphi*v.imphi
      else:            restlo += v.imphi*v.imphi
      if v.implo>=0.0: resthi += v.implo*v.implo
      else:            restlo += v.implo*v.implo
    resthi=  sqrt(resthi)
    restlo= -sqrt(restlo)
    print "Remaining %d impacts after NP %d = %+.6f %+.6f" % (ncons-opt.sumrest, opt.sumrest, resthi, restlo)
    v=vars[isort[opt.sumrest]]
    v.name="#splitline{Remaining systematics}{(quadrature sum)}"
    v.imphi= resthi
    v.implo= restlo
    v.val= -999.0
    v.errlo= v.errhi= 0.0
    ncons=opt.sumrest+1

  bextra= 0.0
  if opt.maxsys>=0:
    ncons=opt.maxsys
    if   opt.haveimp:   restmax= abs(vars[isort[ncons]].impmx)
    elif opt.have_corr: restmax= abs(vars[isort[ncons]].poicor[0])
    else:               restmax= abs(vars[isort[ncons]].poicorfit[0])
    bextra= 1.0
    bmargin += bextra

  if opt.setval:
    m= re.search(r"^(\d+)=(.*)\+([\d.]+)-([\d.]+)$",opt.setval)
    try:
      i,name,resthi,restlo= m.groups()
      i=int(i)
      restlo=float(restlo)
      resthi=float(resthi)
      v=vars[isort[i]]
    except:
      print "Bad spec: --setval",opt.setval
      exit(2)
    v.name=name
    v.imphi=  resthi
    v.implo= -restlo
    v.val= -999.0
    v.errlo= v.errhi= 0.0
    ncons=i+1

  margins= tmargin+bmargin
  cormax *= 1.1
  corscale= pullmax/cormax
  havetxt= bool(opt.txtfile)
  for iplot,ipull in enumerate(xrange(0,ncons,opt.pullsPerPlot)):
    lines= []
    npull= min(opt.pullsPerPlot,ncons-ipull)
    ymax=npull+margins
    jpull= ipull+npull-1
    bottom= 0.08
    if npull<opt.pullsPerPlot: bottom += (1.0-top-bottom) * float(opt.pullsPerPlot-npull)/(opt.pullsPerPlot+margins)
    canvas.SetBottomMargin (bottom)
    zero=     TVectorD(npull)
    posn=     TVectorD(npull)
    widthhi1= TVectorD(npull)
    widthlo1= TVectorD(npull)
    widthhi2= TVectorD(npull)
    widthlo2= TVectorD(npull)
    corveclo1=TVectorD(npull)
    corvechi1=TVectorD(npull)
    corveclo2=TVectorD(npull)
    corvechi2=TVectorD(npull)
    toymeans= TVectorD(npull)
    toyrmslo= TVectorD(npull)
    toyrmshi= TVectorD(npull)
    fitval=   TVectorD(npull)
    fiterrlo= TVectorD(npull)
    fiterrhi= TVectorD(npull)
    pullFrame= ROOT.TH2D ("pullFrame", "", 1, -pullmax, pullmax, npull+2, 0.0, ymax)
    ax= pullFrame.GetYaxis()
    ax.Set (npull+2, array.array("d",[0.0]+[bmargin+i for i in range(npull+1)]+[npull+margins]))
    if   opt.pullsPerPlot>=40: ax.SetLabelSize(0.018)
    elif opt.pullsPerPlot>=15: ax.SetLabelSize(0.027)
    else:                      ax.SetLabelSize(0.040)
    ax.SetTickLength(0.0)
    for i in xrange(npull):
      j= ipull+i
      if opt.haveimp or opt.include_corr==1: j= isort[j]
      v= vars[j]
      if   opt.haveimp:     mucorhi= abs(v.imphi); mucorlo= abs(v.implo)
      elif opt.have_corr:   mucorhi= mucorlo= v.poicor[0]
      elif opt.have_corfit: mucorhi= mucorlo= v.poicorfit[0]
      else:                 mucorhi= mucorlo= 0.0
      ax.SetBinLabel (npull-i+1, niceName(v.name))
      posn[i]= npull-i-0.5+bmargin
      widthhi1[i]= widthlo1[i]= widthhi2[i]= widthlo2[i]= 0.5
      if opt.txtfile:
        corvechi1[i]= abs(corscale*mucorhi)
        corveclo1[i]= abs(corscale*mucorlo)
        if v.havetxt:
          fiterrlo[i]= corscale*v.fiterrlo
          fiterrhi[i]= corscale*v.fiterrhi
      else:
        if opt.haveimp:
          if   not opt.sign_up:
            if   v.imphi*v.implo>0.0:
              if v.imphi>=0.0:
                corvechi1[i]=  v.imphi*corscale
                corvechi2[i]=  v.implo*corscale
                widthlo1[i]= widthhi2[i]= 0.0
              else:
                corveclo1[i]= -v.implo*corscale
                corveclo2[i]= -v.imphi*corscale
                widthhi1[i]= widthlo2[i]= 0.0
            elif v.imphi>=v.implo:
              corvechi1[i]=  max ( v.imphi*corscale,  v.implo*corscale, 0.0)
              corveclo1[i]=  max (-v.imphi*corscale, -v.implo*corscale, 0.0)
              widthlo2[i]= widthhi2[i]= 0.0
            else:
              corvechi2[i]=  max (-v.imphi*corscale, -v.implo*corscale, 0.0)
              corveclo2[i]=  max ( v.imphi*corscale,  v.implo*corscale, 0.0)
              widthlo1[i]= widthhi1[i]= 0.0
          elif not opt.no_sign:
            if v.imphi>=0.0: corvechi1[i]=  v.imphi*corscale
            else:            corvechi2[i]= -v.imphi*corscale
            if v.implo>=0.0: corveclo2[i]=  v.implo*corscale
            else:            corveclo1[i]= -v.implo*corscale
          else:
            corvechi1[i]= corscale*mucorhi
            corveclo1[i]= corscale*mucorlo
        else:
          corvechi1[i]= corscale*mucorhi
          corveclo1[i]= 0.0
          corveclo2[i]= corscale*mucorhi
          corvechi2[i]= 0.0
        fitval[i]=    v.val
        fiterrlo[i]= -v.errlo
        fiterrhi[i]=  v.errhi
      if opt.results_file2:
        toymeans[i]=    v.val2
        toyrmslo[i]= -v.errlo2
        toyrmshi[i]=  v.errhi2
      elif v.ntoys:
        toymeans[i]= v.toymean
        toyrmslo[i]= toyrmshi[i]= v.toyrms
      if i>0 and i<npull and i%nlines==0:
        lines.append (ROOT.TLine (-pullmax, posn[i]+0.5, pullmax, posn[i]+0.5))

    corax= ROOT.TGaxis (-pullmax,ymax,pullmax,ymax,-100.0*cormax,100.0*cormax,510,"-")
    corax.ImportAxisAttributes(pullFrame.GetXaxis())
    # corax.SetTitle("#rho(#mu,P)")

    for lx in [-1.0, 0.0, 1.0]:
      if opt.txtfile and lx!=0.0: continue
      l= ROOT.TLine (lx, 0.0, lx, npull+bmargin)
      l.SetLineWidth(2)
      l.SetLineColor(ROOT.kBlue)
      lines.append(l)
    if verbose>=1:
      print "%-45s %6s %6s %6s %6s" % tuple("NP HI1 LO1 HI2 LO2".split())
      for i in xrange(npull):
        print "%-45s %6.3f %6.3f %6.3f %6.3f" % tuple([vars[isort[ipull+i]].name] + [x[i]/corscale for x in corvechi1, corveclo1, corvechi2, corveclo2])
    corGraph1=ROOT.TGraphAsymmErrors (zero,     posn, corveclo1,corvechi1,widthlo1, widthhi1)
    corGraph2=ROOT.TGraphAsymmErrors (zero,     posn, corveclo2,corvechi2,widthlo2, widthhi2)
    toyGraph= ROOT.TGraphAsymmErrors (toymeans, posn, toyrmslo, toyrmshi, zero,  zero)
    fitGraph= ROOT.TGraphAsymmErrors (fitval,   posn, fiterrlo, fiterrhi, zero,  zero)
    corGraph1.SetFillColor(ROOT.kYellow)
    corGraph1.SetLineColor(0)
    corGraph1.SetLineWidth(0)
    if showcorr==2:
      corGraph2.SetFillStyle(0)
      corGraph2.SetLineColor(ROOT.kYellow)
      corGraph2.SetLineWidth(2)
    else:
      corGraph2.SetFillColor([ROOT.kCyan,ROOT.kYellow][opt.no_sign])
      corGraph2.SetLineColor(0)
      corGraph2.SetLineWidth(0)
    toyGraph.SetMarkerStyle(20)
    toyGraph.SetMarkerSize(1)
    toyGraph.SetLineWidth(2)
    toyGraph.SetLineColor(ROOT.kBlue)
    toyGraph.SetMarkerColor(ROOT.kBlue)
    fitGraph.SetMarkerStyle(24)
    fitGraph.SetMarkerSize((1,0)[havetxt])
    fitGraph.SetLineWidth(2)
    fitGraph.SetLineColor((ROOT.kRed,ROOT.kBlue)[havetxt])
    fitGraph.SetMarkerColor(fitGraph.GetLineColor())
    if opt.txtfile: pullFrame.GetXaxis().SetNdivisions(0)

    leg= ROOT.TLegend (0.5, bottom-(0.08,0.05)[havetxt], 0.97, bottom-(0.04,0.01)[havetxt])
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    if opt.txtfile:
      leg.AddEntry (fitGraph, "Impact on #mu [%]",       "lp")
      leg.AddEntry (corGraph1,"#rho(#mu,P) [%]",         "f")
    elif opt.results_file2 or opt.ntoys:
      leg.AddEntry (fitGraph, "Data fit",            "lp")
      if   opt.results_file2: leg.AddEntry (toyGraph, "Alt fit", "lp")
      elif opt.ntoys:         leg.AddEntry (toyGraph, "Toy mean #pm #sigma", "lp")
    corleg= ROOT.TLegend (left-0.24, 1.02-top-0.02*poilabLines, left-0.04, 1.05-top)
    corleg.SetBorderSize(0)
    corleg.SetLineColor(0)
    corleg.SetFillStyle(0)
    corleg.AddEntry (corGraph1, poilab,      "f")
    if not opt.no_sign:
      corleg.AddEntry (corGraph2, "anti-corr", "f")

    pullleg= ROOT.TLegend (left-0.22, bottom-0.04, left-0.04, bottom)
    pullleg.SetBorderSize(0)
    pullleg.SetLineColor(0)
    pullleg.SetFillStyle(0)
    pullleg.AddEntry (fitGraph, "(#hat{#theta}-#theta_{0})/#sigma_{#theta}", "lp")

    restlab= None
    if opt.maxsys>=0:
      restlab= ROOT.TPaveText(-0.9*pullmax, bmargin-1.0*bextra, 0.9*pullmax, bmargin-0.2*bextra)
      restlab.SetBorderSize(0)
      restlab.SetFillColor(0)
      restlab.SetTextFont(52)
      restlab.AddText("Impact of all other")
      restlab.AddText("nuisance parameters < %.1g%%" % (100.0*restmax))
      lines.append (ROOT.TLine (-pullmax, bmargin, pullmax, bmargin))
    for l in lines: l.SetLineStyle(2)

    if opt.atlasLabel:
      if opt.preliminary or opt.internal:
        AtlasLabel1= ROOT.TPaveText(-0.9*pullmax, ymax-tmargin+0.1, 0.0, ymax-0.35)
      else:
        AtlasLabel1= ROOT.TPaveText(-0.6*pullmax, ymax-tmargin+0.3, 0.0, ymax-0.45)
      AtlasLabel1.SetBorderSize(0)
      AtlasLabel1.SetFillStyle(0)
      AtlasLabel1.SetFillStyle(0)
      AtlasLabel1.SetTextFont(42)
      if   opt.preliminary: AtlasLabel1.AddText("#bf{#it{ATLAS}} Preliminary");
      elif opt.internal:    AtlasLabel1.AddText("#bf{#it{ATLAS}} Internal");
      else:                 AtlasLabel1.AddText("#bf{#it{ATLAS}}");
      AtlasLabel2= ROOT.TPaveText(0.0, ymax-tmargin+0.1, 0.9*pullmax, ymax-0.35)
      AtlasLabel2.SetBorderSize(0)
      AtlasLabel2.SetFillStyle(0)
      AtlasLabel2.SetTextFont(42)
      AtlasLabel2.SetTextAlign(12)
      if opt.e7TeV: AtlasLabel2.AddText("#sqrt{s} = 7 TeV #scale[0.6]{#int}Ldt = 4.5 fb^{-1}");
      if opt.e8TeV: AtlasLabel2.AddText("#sqrt{s} = 8 TeV #scale[0.6]{#int}Ldt = 20.3 fb^{-1}");

    pullFrame.Draw()
    if showcorr:
      corax.Draw()
      corGraph1.Draw("2")
      if not opt.no_sign or not opt.sign_up or showcorr==2: corGraph2.Draw("2")
    if opt.ntoys or opt.results_file2: toyGraph.Draw("ep")
    fitGraph.Draw("ep")
    for o in lines: o.Draw()
    leg.Draw()
    if showcorr and not opt.txtfile: corleg.Draw()
    if not opt.txtfile: pullleg.Draw()
    if restlab: restlab.Draw()
    if opt.atlasLabel:
      AtlasLabel1.Draw();
      AtlasLabel2.Draw();
    PrintPlot(canvas,pullFile)

    if opt.no_latex: continue

    draw= []
    name= "pulls_"+str(iplot)
    eps= ROOT.TCanvas (name, name, 700, 500)
    eps.SetFillStyle(0);
    top=0.05
    bottom= 0.08
    if npull<opt.pullsPerPlot: bottom += (1.0-top-bottom) * float(opt.pullsPerPlot-npull)/(opt.pullsPerPlot+2)
    eps.SetLeftMargin  (0.5)
    eps.SetRightMargin (0.005)
    eps.SetTopMargin   (top)
    eps.SetBottomMargin(bottom)
    if opt.pullsPerPlot>=40: ax.SetLabelSize(0.03)
    else:                    ax.SetLabelSize(0.04)
    ax.SetLabelSize(0.03)
    corax.SetLabelSize(0.03)
    leg.SetX1NDC(0.6)
    leg.SetX2NDC(0.995)
    corleg.SetX1NDC(0.27)
    corleg.SetX2NDC(0.47)
    corleg.SetY1NDC(0.995-top)
    corleg.SetY2NDC(1.045-top)
    pullleg.SetX1NDC(0.27)
    pullleg.SetX2NDC(0.47)
    draw.append(pullFrame)
    if showcorr:
      draw.append(corax)
      draw.append([corGraph1,"2"])
      if not opt.no_sign or not opt.sign_up or showcorr==2: draw.append([corGraph2,"2"])
    if opt.ntoys or opt.results_file2: draw.append([toyGraph,"ep"])
    draw.append([fitGraph,"ep"])
    draw.extend(lines)
    draw.append(leg)
    if showcorr and not opt.txtfile: draw.append(corleg)
    if not opt.txtfile: draw.append(pullleg)
    PlotEps (name, draw, eps)

    lab=""
    if iplot==0:       lab += "\n\\label{fig:pulls}"
    if jpull>=ncons-1: lab += "\n\\label{fig:pulls_end}"

    print >>ptex, ptex_template % ("figures/%s/%s" % (opt.outname,name), lab)

  PlotFile(canvas,pullFile,True)
  if not opt.no_latex: del ptex


def RunLaTeX (texfile_main,
              texfile_include1, texfile_include2, texfile_include3, texfile_include4,
              texfile1, texfile2, texfile3, texfile4,
              figdir, fig_path, pdffile):

  dvifile= re.sub(r"\.tex$","",texfile_main) + ".dvi"
  parent_dir= gSystem.DirName(fig_path)
  cc= ("","#")
  comm= cc[texfile1==texfile_include1]
  runcmd ("""\
set -x
%sln -nfs '%s' '%s'
%sln -nfs '%s' '%s'
%sln -nfs '%s' '%s'
%sln -nfs '%s' '%s'
%sln -nfs '%s' '%s'
%sln -nfs '%s' '%s'
mkdir -p '%s'
ln -nfs '%s' '%s'
latex '%s' && latex '%s' && latex '%s' && latex '%s' && dvipdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress '%s' '%s'
rm '%s'
rmdir '%s'
""" % (comm, texfile1%"unconstrained", texfile_include1%"unconstrained",
       comm, texfile1%"constrained",   texfile_include1%"constrained",
       comm, texfile1%"gamma",         texfile_include1%"gamma",
       cc[texfile2==texfile_include2], texfile2, texfile_include2,
       cc[texfile3==texfile_include3], texfile3, texfile_include3,
       cc[texfile4==texfile_include4], texfile4, texfile_include4,
       parent_dir, figdir, fig_path,
       texfile_main, texfile_main, texfile_main, texfile_main,
       dvifile, pdffile, fig_path, parent_dir))


def MergeTree (f, pll, name="pll"):
  sadd= f.Get (name)
  if not sadd:
    print "No object %s in file %s" % (name, f.GetName())
    return None
  print "Read '%s' from file %s" % (name, f.GetName())
  if pll==None:
    pll= sadd
  else:
    tadd= ROOT.TList(sadd)
    pll.Merge (tadd)
  return pll


def varInfo (v, isres2=False):
  name= v.GetName()
  val=  v.getVal()
  err= -1.0
  errlo=0.0
  errhi=0.0
  errtype= "C"
  if   v.hasAsymError(0): errqual=2
  elif v.hasAsymError():  errqual=1
  else:                   errqual=0
  if not v.isConstant():
    if       v.hasAsymError():
      errlo= v.getAsymErrorLo()
      errhi= v.getAsymErrorHi()
      if  errlo>=0.0:
        err =  errhi
        errtype= "+"
      elif errhi<=0.0:
        err = -errlo
        errtype= "-"
      else:
        err = 0.5 * (errhi - errlo)
        errtype= " "
    elif v.hasError():
      err = v.getError()
      errtype= "*"
  if opt.use_migrad_errs:
    if errlo>=0: errlo= -err
    if errhi<=0: errhi=  err

  symm= (err<=0.0 or abs ((errhi+errlo)/err) < 0.05)
  if err>0.0: nf= max (0, opt.sf-1-int(math.floor(math.log10(err*1.05))))
  else:       nf= 2

  v2= ROOT.RooRealVar (name, name, -999.0)
  if isres2:
    return Struct (name=name, have_data=False, havetxt=False, haveimp=0, havecont=0, haveval=False, haveval2=True,
                   var=v2, val=-999.0, err=-1.0, errlo=0.0, errhi=0.0, errqual=-1, errtype="X", symm=False, nf=0,
                   var2=v, val2=val, err2=err, errlo2=errlo, errhi2=errhi, errqual2=errqual, errtype2=errtype, symm2=symm, nf2=nf,
                   imphi=0.0, implo=0.0, impflag="")
  else:
    return Struct (name=name, have_data=False, havetxt=False, haveimp=0, havecont=0, haveval=True, haveval2=False,
                   var=v, val=val, err=err, errlo=errlo, errhi=errhi, errqual=errqual, errtype=errtype, symm=symm, nf=nf,
                   var2=v2, val2=-999.0, err2=-1.0, errlo2=0.0, errhi2=0.0, errqual2=-1, errtype2="X", symm2=False, nf2=0,
                   imphi=0.0, implo=0.0, impflag="")


def PrintPlot(canvas,plotFile):
  oldErrorIgnoreLevel= ROOT.gErrorIgnoreLevel
  if verbose<1: ROOT.gErrorIgnoreLevel= ROOT.kWarning
  canvas.Print (plotFile)
  ROOT.gErrorIgnoreLevel= oldErrorIgnoreLevel
  if interactive:
    # pause until user double-clicks or presses a key in the window before going to next plot
    canvas.Update()
    canvas.WaitPrimitive()
  canvas.Clear()


def PlotFile(canvas,plotFile,end=False):
  oldErrorIgnoreLevel= ROOT.gErrorIgnoreLevel
  if verbose<1 and not end: ROOT.gErrorIgnoreLevel= ROOT.kWarning
  if not end:
    canvas.Print (plotFile+"[")
  else:
    canvas.Print (plotFile+"]")
    canvas.IsA().Destructor(canvas)
  ROOT.gErrorIgnoreLevel= oldErrorIgnoreLevel


def PlotPdf (canvas,plotFile,page,nx=2,ny=2):
  if nx*ny>1: canvas.Divide (nx, ny, 0.001, 0.001)
  for i,draw in enumerate(page):
    canvas.cd(i+1)
    gPad.SetLogy()
    for dl in draw:
      if isinstance(dl,list): d,o= dl
      else:                   d,o= [dl,""]
      d.Draw(o)
  canvas.cd()
  oldErrorIgnoreLevel= ROOT.gErrorIgnoreLevel
  if verbose<1: ROOT.gErrorIgnoreLevel= ROOT.kWarning
  canvas.Print (plotFile)
  ROOT.gErrorIgnoreLevel= oldErrorIgnoreLevel
  if interactive:
    # pause until user double-clicks or presses a key in the window before going to next plot
    canvas.Update()
    canvas.WaitPrimitive()
  canvas.Clear()


def PlotEps (name, draw, use_canvas=None):
  global tmpdir
  if not tmpdir:
    tmp= gSystem.TempDirectory()
    try:
      tmpdir= tempfile.mkdtemp ("","npnote_",tmp)
    except:
      print "Could not create temporary directory in",tmp
      return
    if verbose>=0: print "Write EPS files to",tmpdir

  oldPad= gPad.GetPad(0)

  eps= use_canvas
  ROOT.gROOT.SetBatch()
  if not eps:
    eps= ROOT.TCanvas (name, name, 700, 700)
    eps.SetFillStyle(0);
    eps.SetLogy()
    eps.SetTopMargin(0.07)
    eps.SetBottomMargin(0.04)
    eps.SetLeftMargin(0.06)
    eps.SetRightMargin(0.01)
  eps.cd()

  for dl in draw:
    if isinstance(dl,list): d,o= dl
    else:                   d,o= [dl,""]
    d.Draw(o)
  if interactive: ROOT.gROOT.SetBatch(0)

  oldErrorIgnoreLevel= ROOT.gErrorIgnoreLevel
  if verbose<1: ROOT.gErrorIgnoreLevel= ROOT.kWarning
  eps.Print (tmpdir+"/"+name+".eps")
  ROOT.gErrorIgnoreLevel= oldErrorIgnoreLevel
  eps.Clear()
  eps.IsA().Destructor(eps)
  if oldPad: oldPad.cd()


def PlotNP (nt, varname, name, dval, errlo, errhi, cut, maxEnt, allEntries):
  if not nt: return 0, None, None
  haveWeight=nt.GetBranch("weight")
  haveAlt=nt.GetBranch("is_alt")
  if allEntries:
    npcut= cut
  else:
    npcut= "!is_data"
    if haveAlt:               npcut += " && !is_alt"
    if opt.plotMuHatPositive: npcut += " && mu>=0.0"
    if cut!="":               npcut += " && (%s)" % cut
    if haveWeight:            npcut  = "(%s)*weight" % npcut
  if verbose>=1: print "Ntuple var",name,"selector:",npcut

  draw= []
  if maxEnt<0: maxEnt= 1000000000
  nv= nt.Draw (varname, npcut, "goff", maxEnt)
# if varname=="mu": nv= nt.Draw (varname+">>(100,0,5)", npcut, "goff", maxEnt)
  hn= nt.GetHistogram()
  if nv==0 or not hn: return 0, None, None
  hn= hn.Clone(name)
  hn.SetTitle(name)
  hn.GetXaxis().SetTitle()
  draw.append(hn)
  xnlo=hn.GetXaxis().GetXmin()
  xnhi=hn.GetXaxis().GetXmax()
  xnb=(xnhi-xnlo)/hn.GetNbinsX()
  sumw= hn.GetSum()
  if verbose>=1: print "%s: sumw=%g, mu=%g, err=(%g,%g), range=%g:%g" % (varname, sumw, dval, errlo, errhi, xnlo, xnhi)
  if errlo<0.0 and errhi>0.0:
    c= sumw*xnb* (0.5/errhi - 0.5/errlo)/sqrt(2.0*pi)
    gaus= ROOT.TF1 ("gaus", "[0]*exp(-0.5*((x-[1])/(x<[1]?[2]:[3]))**2)", xnlo, xnhi)
    gaus.SetParNames     ("Constant", "Data fit", "Minos -err", "Minos +err")
    gaus.SetParameters   (c,          dval,       -errlo,       errhi)
    gaus.FixParameter (0, c)
    gaus.SetNpx(3*hn.GetNbinsX())
    gaus.SetLineColor(ROOT.kRed)
    draw.append([gaus,"lsame"])
    if hn.GetListOfFunctions(): hn.GetListOfFunctions().Add(gaus)
  elif dval>=xnlo and dval<=xnhi:
    dline= ROOT.TLine (dval, hn.GetMinimum(), dval, 0.8*hn.GetMaximum())
    dline.SetLineWidth(3)
    dline.SetLineColor(ROOT.kRed)
    draw.append(dline)
  return nv, hn, draw


def DrawLatex (x, y, text, ndc=True, align=31, size=0.03, color=None):
  if color==None: color= ROOT.kBlack
  t= ROOT.TLatex(x,y,text)
  t.SetTextAlign(align)
  t.SetTextSize(size)
  if ndc: t.SetNDC()
  t.SetTextColor(color)
  return t


# Parse arguments and run main program
exit (process (parseArgs()))
