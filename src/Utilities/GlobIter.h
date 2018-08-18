//=====================================================================-*-C++-*-
// File: $Id: GlobIter.h 727260 2016-03-02 15:10:46Z adye $
//==============================================================================

#ifndef GlobIter_h
#define GlobIter_h

#include <iostream>
#include <vector>
#include <algorithm>

#include "TError.h"
#include "TSystem.h"
#include "TString.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "TClass.h"
#include "TKey.h"
#include "TFile.h"

//==============================================================================
// FileGlobIter
//==============================================================================

namespace GlobIter {
  enum ErrorMode { kMessage, kIgnore, kReturn };
  int debug = 0;
}

class FileGlobIter {

private:
  std::vector<std::string> _files;
  TString _spec, _pre, _suf, _current;
  GlobIter::ErrorMode _errorMode;
  int _depth, _nret;
  size_t _ifile;
  bool _isname, _started;
  FileGlobIter* _subIter;
  TFile* _file;

  void Setup() {
    _started = false;
    _subIter = 0;
    _ifile = 0;
    _nret = 0;
    _file= 0;

    if (GlobIter::debug>=1) std::cout << "FileGlobIter " << _depth << ": " << _spec << std::endl;
    if (_depth>50) return;  // just in case
    TString path = _spec;
    if (_depth==0) gSystem->ExpandPathName(path);
    static const char *specials = "[]*?";
    int wild = path.First(specials);
    if (wild == kNPOS) {
      if (GlobIter::debug>=1) std::cout << "FileGlobIter " << _depth << ": no wild " << path << std::endl;
      _isname= true;
      if (gSystem->AccessPathName(path)) return;  // don't save anything if a non-wildcarded file doesn't exist
      _files.push_back(path.Data());
      return;
    }

    int next= -1, slash= 0;
    do {
      slash= next;
      next = path.Index('/',slash+1);
    } while (next != kNPOS && next < wild);
    _isname = (next == kNPOS);

    TString dir, sub;
    if (slash>0) {
      dir = path(0,slash);
      _pre= dir;
      _pre += '/';
    } else
      dir = ".";

    if (next != kNPOS) {
      _suf= path(next,path.Length()-next);
      sub = path(slash+1,next-slash-1);
    } else {  // file name is wildcarded
      sub = path(slash+1,path.Length()-slash-1);
    }

    if (GlobIter::debug>=1) std::cout << "FileGlobIter " << _depth << ": wild " << dir << " " << sub << std::endl;

    void* dirlist = gSystem->OpenDirectory (dir);
    if (!dirlist) return;
    // Create regular expression matcher using wildcard expression.
    // Note that TRegexp does not escape ^, $, or + so those characters should be escaped by the user.
    TRegexp re(sub,kTRUE);
    bool sel= (sub!="*");
    const char* file;
    while ((file = gSystem->GetDirEntry(dirlist))) {
      if (!strcmp(file,".") || !strcmp(file,"..")) continue;
      if (sel) {
        TString s = file;
        if (s != sub && s.Index(re) == kNPOS) continue;
      }
      _files.push_back (file);
    }
    gSystem->FreeDirectory(dirlist);
#ifndef __CINT__
    // sort doesn't seem to work from CINT, so we return unsorted
    std::sort (_files.begin(), _files.end());
#endif
  }

  const char* UseFile() {
    _nret++;
    return _current.Data();
  }

  const char* UseFile (const TString& file) {
    _current= file;
    return UseFile();
  }

  FileGlobIter (const char* spec, int depth)
    : _spec(spec), _errorMode(GlobIter::kIgnore), _depth(depth) { Setup(); }

public:

  FileGlobIter (const char* spec, GlobIter::ErrorMode errorMode = GlobIter::kMessage)
    : _spec(spec), _errorMode(errorMode),         _depth(0)     { Setup(); }

  const char* Next() {
    delete _file; _file= 0;
    _started= true;
    for (;;) {
      if (_subIter) {
        const char* f = _subIter->Next();
        if (f) return UseFile(f);
        delete _subIter;
        _subIter= 0;
      }
      if (_ifile >= _files.size()) {
        if (_nret == 0) {  // no files matched
          if      (_errorMode == GlobIter::kMessage) Error ("FileGlobIter", "File %s not found", _spec.Data());
          else if (_errorMode == GlobIter::kReturn)  return UseFile(_spec);
        }
        _started= false;
        return 0;
      }
      _current = _pre + _files[_ifile++].c_str() + _suf;
      if (_isname) return UseFile();
      _subIter = new FileGlobIter(_current, _depth+1);
    }
  }

  const char* operator() ()     { return Next(); }
  const char* operator*() const { return _started ? _current.Data() : 0; }

  void Reset() {
    delete _subIter; _subIter= 0;
    delete _file; _file= 0;
    _started= false;
    _files.clear();
    _pre.Clear();
    _suf.Clear();
    Setup();
  }

  const char* GetSpec() const { return _spec.Data(); }

  int NumFiles() const { return _nret; }  // number of files returned so far

  TFile* GetFile (Option_t* option="") {
    if (_file) return _file;
    _file= TFile::Open (_current.Data(), option);
    return _file;
  }

  TFile* NextFile (Option_t* option="") { return Next() ? GetFile(option) : 0; }

  ~FileGlobIter() {
    delete _subIter;
    delete _file;
  }
};

//==============================================================================
// TDirectoryGlobIter
//==============================================================================

class TDirectoryGlobIter {

private:
  TDirectory* _dir;
  TString _spec, _dirname, _pre, _suf, _current;
  GlobIter::ErrorMode _errorMode;
  int _depth, _nret;
  bool _isname, _started, _single;
  TKey* _key;
  TDirectoryGlobIter* _subIter;
  TIter* _iter;
  TRegexp* _re;
  TObject* _obj;

  void Setup() {
    _started = false;
    _single = false;
    _key = 0;
    _subIter = 0;
    _iter = 0;
    _re = 0;
    _nret = 0;
    _obj = 0;

    if (GlobIter::debug>=1) std::cout << "TdirectoryGlobIter " << _depth << ": " << _dir->GetPath() << " " << _spec << std::endl;

    if (_depth>50) return;  // just in case
    TString path = _spec;
    static const char *specials = "[]*?";
    int wild = path.First(specials);

    TString dir, sub;
    TDirectory* subdir;

    if (wild == kNPOS) {
      if (GlobIter::debug>=1) std::cout << "TdirectoryGlobIter " << _depth << ": no wild " << _dir->GetPath() << " " << path << std::endl;
      _isname= true;
      int slash= path.Last('/');
      if (slash != kNPOS) {
        dir = path(0,slash);
        _pre= dir;
        sub = path(slash+1,path.Length());
        subdir = _dir->GetDirectory(dir);
      } else {
        subdir = _dir;
        sub = path;
      }
      if (subdir) _key = subdir->FindKey(sub);
      _single= _key;
      return;
    }

    int next= -1, slash= 0;
    do {
      slash= next;
      next = path.Index('/',slash+1);
    } while (next != kNPOS && next < wild);
    _isname = (next == kNPOS);

    if (slash>0) {
      dir = path(0,slash);
      _pre= dir;
      subdir = _dir->GetDirectory (dir);
      if (!subdir) return;
    } else
      subdir = _dir;

    if (next != kNPOS) {
      _suf= path(next+1,path.Length()-next-1);
      sub = path(slash+1,next-slash-1);
    } else {  // file name is wildcarded
      sub = path(slash+1,path.Length()-slash-1);
    }

    if (GlobIter::debug>=1) std::cout << "TdirectoryGlobIter " << _depth << ": wild " << subdir->GetPath() << " " << sub << std::endl;

    TList* dirl= subdir->GetListOfKeys();
    if (!dirl) return;
    _iter = new TIter(dirl);
    if (sub!="*") _re = new TRegexp(sub,kTRUE);
    else          _re= 0;
  }

  const char* UseKey (TKey* key, const char* name=0) {
    _key= key;
    _current = _dirname;
    if (!_current.IsNull()) _current += "/";
    if (!_pre.IsNull()) {
      _current += _pre;
      _current += "/";
    }
    _current += name ? name : key->GetName();
    _nret++;
    return _current.Data();
  }

  TDirectoryGlobIter (TDirectory* dir, const char* spec, const char* dirname, int depth)
    : _dir(dir), _spec(spec), _dirname(dirname), _errorMode(GlobIter::kIgnore), _depth(depth) { Setup(); }

public:

  TDirectoryGlobIter (TDirectory* dir, const char* spec, GlobIter::ErrorMode errorMode = GlobIter::kMessage)
    : _dir(dir), _spec(spec), _errorMode(errorMode), _depth(0) { Setup(); }

  const char* Next() {
    delete _obj; _obj= 0;
    _started= true;
    for (;;) {
      if (_subIter) {
        const char* f = _subIter->Next();
        if (f) return UseKey (_subIter->GetKey(), f);
        delete _subIter;
        _subIter= 0;
      }
      if (_iter) {
        TKey* k;
        while ((k= dynamic_cast<TKey*>(_iter->Next()))) {
          if (_re && TString(k->GetName()).Index(*_re) == kNPOS) continue;
          const TKey* ktop= k->GetMotherDir()->GetKey (k->GetName());  // check this is the highest cycle
          if (ktop && k->GetCycle() != ktop->GetCycle()) continue;     // (would be more efficient to compare with last name)
          if (_isname) return UseKey (k);
          const TClass* c= TClass::GetClass (k->GetClassName());
          if (!(c && c->InheritsFrom (TDirectory::Class()))) continue;
          TObject* obj= k->ReadObj();
          TDirectory* sub= dynamic_cast<TDirectory*>(obj);
          if (!sub) {
            delete obj;
            continue;
          }
          _subIter = new TDirectoryGlobIter (sub, _suf, k->GetName(), _depth+1);
          break;
        }
        if (_subIter) continue;
      }
      if (_iter || !_single) {
        if (_nret == 0) {  // no keys matched
          if      (_errorMode == GlobIter::kMessage)
            Error ("FileGlobIter", "Key %s not found in %s", _spec.Data(), _dir->GetPath());
          else if (_errorMode == GlobIter::kReturn) {
            _current = _spec;
            return _current.Data();
          }
        }
        _started= false;
        return 0;
      }
      _single= false;
      return UseKey (_key);
    }
  }

  const char* operator() ()     { return Next(); }
  TKey*       NextKey()         { return Next()   ? _key            : 0; }

  const char* operator*() const { return _started ? _current.Data() : 0; }
  TKey*       GetKey()    const { return _started ? _key            : 0; }
  const char* GetName()   const { const TKey* k= GetKey(); return k ? k->GetName() : 0; }

  void Reset() {
    delete _subIter; _subIter= 0;
    delete _iter;    _iter= 0;
    delete _re;      _re= 0;
    delete _obj;     _obj= 0;
    _started= false;
    _key= 0;
    _pre.Clear();
    _suf.Clear();
    Setup();
  }

  const char* GetSpec() const { return _spec.Data(); }
  TDirectory* GetBaseDirectory() const { return _dir; }

  TDirectory* GetDirectory() const {
    TKey* k= GetKey();
    if (!k) return 0;
    return k->GetMotherDir();
  }

  const char* GetPath() const {
    TKey* k= GetKey();
    if (!k) return 0;
    TDirectory* dir= k->GetMotherDir();
    if (!dir) return 0;
    static TString d;
    d= dir->GetPath();
    if (d.IsNull()) return 0;
    if (!d.EndsWith("/")) d += "/";
    d += k->GetName();
    return d.Data();
  }

  int NumKeys() const { return _nret; }  // number of keys returned so far

  TObject* GetObj  (const TClass* requiredClass=0) {
    if (_obj) return _obj;
    if (requiredClass) {
      const TClass* c= TClass::GetClass (_key->GetClassName());
      if (!(c && c->InheritsFrom (requiredClass))) return 0;
    }
    _obj= _key->ReadObj();
    return _obj;
  }

  TObject* NextObj (const TClass* requiredClass=0) { return Next() ? GetObj(requiredClass) : 0; }

  ~TDirectoryGlobIter() {
    delete _subIter;
    delete _iter;
    delete _re;
    delete _obj;
    if (_depth>0) delete _dir;
  }
};

#endif
