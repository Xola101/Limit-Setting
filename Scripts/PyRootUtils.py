# $Id: PyRootUtils.py 676815 2015-06-19 18:56:19Z adye $
# Various PyROOT utilities
# Importing this module loads RooFit without any startup message.
# It also makes RooAbsCollection iterable.

import sys, os
import ROOT
__all__ = ["Struct", "cond", "isnan", "copysign", "runcmd", "CompileCpp", "FileGlobIter", "TDirectoryGlobIter"]


doneInitRooFit=False
def initRoot():
  """Inhibit annoying RooFit startup message"""
  global doneInitRooFit
  if doneInitRooFit: return
  doneInitRooFit= True

  saveErrorIgnoreLevel= ROOT.gErrorIgnoreLevel
  ROOT.gErrorIgnoreLevel= ROOT.kFatal     # Various ROOT errors from rdbuf(), but it looks like we can ignore them
  devnull= ROOT.std.ofstream (os.devnull) # Open /dev/null
  coutbuf= ROOT.std.cout.rdbuf()          # save old buf
  ROOT.std.cout.rdbuf (devnull.rdbuf())   # redirect cout to /dev/null
  ROOT.gErrorIgnoreLevel= saveErrorIgnoreLevel
  ROOT.RooFit                             # Load RooFit
  ROOT.std.cout.rdbuf (coutbuf)           # restore old cout

# PyROOT in ROOT 6 doesn't print the RooFit startup message
# In fact, does any version of PyROOT? The following gave an error in ROOT 5.34.19+, so don't do it.
#if ROOT.gROOT.GetVersionCode() < 393216: initRoot()


class Struct(dict):
  """
  C struct-like class. Usage:
      obj=Struct(a=1,b=2)
      obj.b
      obj.c=3
  See http://code.activestate.com/recipes/52308/#c4
  """
  def __init__(self,*a,**kw):
    dict.__init__(self,*a,**kw)
    self.__dict__= self
  def __add__(self,other):
    new=Struct(self)
    new.update(other)
    return new

def cond(test,iftrue,iffalse):
  """
  Equivalent to "IFTRUE if TEST else IFFALSE" without shortcutting,
  but works before Python 2.5 ("cond" was a disfavoured proposal for PEP-308)
  """
  if test: return iftrue
  else:    return iffalse

try:
  from math import isnan
except ImportError:
  def isnan(x):
    """Check if the float x is a NaN (not a number). Like math.isnan(x), but also works before Python 2.6"""
    return x != x

try:
  from math import copysign
except ImportError:
  def copysign(x,y):
    """Return x with the sign of y. like math.copysign(x,y), but also works before Python 2.6"""
    y= float(y)
    if y<0.0:          return -abs(float(x))
    if y>0.0:          return  abs(float(x))
    if str(y)[0]=='-': return -abs(float(x))   # -0.0 returns x<0
    else:              return  abs(float(x))

def runcmd (cmd):
  sys.stdout.flush()
  stat= os.system (cmd)
  if stat==0: return 0
  exstat= (stat>>8)&255
  if exstat==0: exstat= (stat&255)|128
  return exstat


# Enable Python-style iteration over RooAbsCollection
def TIterator_next(self):
  v=self.Next()
  if v: return v
  raise StopIteration()
ROOT.TIterator.next= TIterator_next  # handles all TIterators, including RooLinkedListIter, returned by RooAbsCollection.iterator()
ROOT.TIterator.__iter__= lambda self: self
ROOT.RooAbsCollection.__iter__= lambda self: self.iterator()


doneCompileCpp= {}
def CompileCpp (filename, reload=False, opt="k", library_name="", build_dir="", dirmode=0):
  if not reload and filename in doneCompileCpp: return True
  filespec= ROOT.gSystem.Which (ROOT.TROOT.GetMacroPath(), filename, ROOT.kReadPermission)
  if filespec=="": filespec= filename   # get error message from CompileMacro
  if not ROOT.gSystem.CompileMacro (filespec, opt, library_name, build_dir, dirmode): return False
  doneCompileCpp[filename]= True
  return True


def FileGlobIter(*args,**kw):
  """File glob iterator. Probably it's better to use glob.glob."""
  if not CompileCpp("GlobIter.h"): return None
  ROOT.FileGlobIter.__iter__= lambda self: self
  ROOT.FileGlobIter.next= TIterator_next
  return ROOT.FileGlobIter(*args,**kw)


def TDirectoryGlobIter(*args,**kw):
  """TDirectory glob iterator."""
  if not CompileCpp("GlobIter.h"): return None
  ROOT.TDirectoryGlobIter.__iter__= lambda self: self
  ROOT.TDirectoryGlobIter.next= TIterator_next
  return ROOT.TDirectoryGlobIter(*args,**kw)
