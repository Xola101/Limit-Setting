# $Id: RooPyUtils.py 578388 2014-01-14 18:36:40Z adye $
# Interface to C++ routines in RooPyUtils.cc.

import ROOT

doneCompileCpp= {}
def CompileCpp (filename, reload=False, opt="k", library_name="", build_dir="", dirmode=0):
  if not reload and filename in doneCompileCpp: return True
  if not ROOT.gSystem.CompileMacro (filename, opt, library_name, build_dir, dirmode): return False
  doneCompileCpp[filename]= True
  return True


def Corr (nt, vars, cut="", maxEnt=-1, verbose=0):
  """
  Return correlation matrix between variables, vars, in ntuple, nt.
  This implementation calls a C++ routine compiled by ACLiC.
  Time to scan 200 variables in 1987 entries: 8 sec.
  The fastest way without using ACLiC takes 144 sec.
  """
  if not CompileCpp("RooPyUtils.cc"): return None
  nvars= len(vars)
  cor= ROOT.TMatrixD (nvars,nvars)
  namevec= ROOT.vector("string")(nvars)
  for i in xrange(nvars): namevec[i]= vars[i]
  ok= ROOT.RooPyUtils.Corr (nt, namevec, cor, cut, maxEnt, verbose)
  if not ok: return None
  return [[cor(i,j) for j in xrange(nvars)] for i in xrange(nvars)]
