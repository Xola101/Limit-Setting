// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__WorkspaceCalculator

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h"

// Header files passed via #pragma extra_include

namespace RooStats {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *RooStats_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("RooStats", 0 /*version*/, "", 42,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &RooStats_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *RooStats_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static void *new_TOwnedList(void *p = 0);
   static void *newArray_TOwnedList(Long_t size, void *p);
   static void delete_TOwnedList(void *p);
   static void deleteArray_TOwnedList(void *p);
   static void destruct_TOwnedList(void *p);
   static Long64_t merge_TOwnedList(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TOwnedList*)
   {
      ::TOwnedList *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TOwnedList >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TOwnedList", ::TOwnedList::Class_Version(), "", 63,
                  typeid(::TOwnedList), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TOwnedList::Dictionary, isa_proxy, 4,
                  sizeof(::TOwnedList) );
      instance.SetNew(&new_TOwnedList);
      instance.SetNewArray(&newArray_TOwnedList);
      instance.SetDelete(&delete_TOwnedList);
      instance.SetDeleteArray(&deleteArray_TOwnedList);
      instance.SetDestructor(&destruct_TOwnedList);
      instance.SetMerge(&merge_TOwnedList);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TOwnedList*)
   {
      return GenerateInitInstanceLocal((::TOwnedList*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TOwnedList*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WorkspaceCalculator(void *p);
   static void deleteArray_WorkspaceCalculator(void *p);
   static void destruct_WorkspaceCalculator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WorkspaceCalculator*)
   {
      ::WorkspaceCalculator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WorkspaceCalculator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WorkspaceCalculator", ::WorkspaceCalculator::Class_Version(), "", 75,
                  typeid(::WorkspaceCalculator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WorkspaceCalculator::Dictionary, isa_proxy, 4,
                  sizeof(::WorkspaceCalculator) );
      instance.SetDelete(&delete_WorkspaceCalculator);
      instance.SetDeleteArray(&deleteArray_WorkspaceCalculator);
      instance.SetDestructor(&destruct_WorkspaceCalculator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WorkspaceCalculator*)
   {
      return GenerateInitInstanceLocal((::WorkspaceCalculator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WorkspaceCalculator*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TOwnedList::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TOwnedList::Class_Name()
{
   return "TOwnedList";
}

//______________________________________________________________________________
const char *TOwnedList::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOwnedList*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TOwnedList::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOwnedList*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TOwnedList::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOwnedList*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TOwnedList::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOwnedList*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WorkspaceCalculator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WorkspaceCalculator::Class_Name()
{
   return "WorkspaceCalculator";
}

//______________________________________________________________________________
const char *WorkspaceCalculator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WorkspaceCalculator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WorkspaceCalculator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WorkspaceCalculator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WorkspaceCalculator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WorkspaceCalculator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WorkspaceCalculator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WorkspaceCalculator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TOwnedList::Streamer(TBuffer &R__b)
{
   // Stream an object of class TOwnedList.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TOwnedList::Class(),this);
   } else {
      R__b.WriteClassBuffer(TOwnedList::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TOwnedList(void *p) {
      return  p ? new(p) ::TOwnedList : new ::TOwnedList;
   }
   static void *newArray_TOwnedList(Long_t nElements, void *p) {
      return p ? new(p) ::TOwnedList[nElements] : new ::TOwnedList[nElements];
   }
   // Wrapper around operator delete
   static void delete_TOwnedList(void *p) {
      delete ((::TOwnedList*)p);
   }
   static void deleteArray_TOwnedList(void *p) {
      delete [] ((::TOwnedList*)p);
   }
   static void destruct_TOwnedList(void *p) {
      typedef ::TOwnedList current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TOwnedList(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TOwnedList*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TOwnedList

//______________________________________________________________________________
void WorkspaceCalculator::Streamer(TBuffer &R__b)
{
   // Stream an object of class WorkspaceCalculator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WorkspaceCalculator::Class(),this);
   } else {
      R__b.WriteClassBuffer(WorkspaceCalculator::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WorkspaceCalculator(void *p) {
      delete ((::WorkspaceCalculator*)p);
   }
   static void deleteArray_WorkspaceCalculator(void *p) {
      delete [] ((::WorkspaceCalculator*)p);
   }
   static void destruct_WorkspaceCalculator(void *p) {
      typedef ::WorkspaceCalculator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WorkspaceCalculator

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 210,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_libWorkspaceCalculator_Impl() {
    static const char* headers[] = {
"/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/local/include",
"/home/nicolin/CERN/LimitSetting/bin/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libWorkspaceCalculator dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate(R"ATTRDUMP(file_name@@@/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h")))  TOwnedList;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(implements the workspace calculator)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h")))  WorkspaceCalculator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libWorkspaceCalculator dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/nicolin/CERN/LimitSetting/src/WorkspaceCalculator.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"NaN", payloadCode, "@",
"Strtod", payloadCode, "@",
"Strtol", payloadCode, "@",
"TOwnedList", payloadCode, "@",
"TOwnedList::fgIsA", payloadCode, "@",
"WorkspaceCalculator", payloadCode, "@",
"WorkspaceCalculator::OptimizeBits", payloadCode, "@",
"WorkspaceCalculator::default_optimize", payloadCode, "@",
"WorkspaceCalculator::fgIsA", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libWorkspaceCalculator",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libWorkspaceCalculator_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libWorkspaceCalculator_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libWorkspaceCalculator() {
  TriggerDictionaryInitialization_libWorkspaceCalculator_Impl();
}
