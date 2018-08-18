// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__FitCalculator

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
#include "/home/nicolin/CERN/LimitSetting/src/FitCalculator.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_FitCalculator(void *p = 0);
   static void *newArray_FitCalculator(Long_t size, void *p);
   static void delete_FitCalculator(void *p);
   static void deleteArray_FitCalculator(void *p);
   static void destruct_FitCalculator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FitCalculator*)
   {
      ::FitCalculator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FitCalculator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FitCalculator", ::FitCalculator::Class_Version(), "", 14,
                  typeid(::FitCalculator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FitCalculator::Dictionary, isa_proxy, 4,
                  sizeof(::FitCalculator) );
      instance.SetNew(&new_FitCalculator);
      instance.SetNewArray(&newArray_FitCalculator);
      instance.SetDelete(&delete_FitCalculator);
      instance.SetDeleteArray(&deleteArray_FitCalculator);
      instance.SetDestructor(&destruct_FitCalculator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FitCalculator*)
   {
      return GenerateInitInstanceLocal((::FitCalculator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FitCalculator*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr FitCalculator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FitCalculator::Class_Name()
{
   return "FitCalculator";
}

//______________________________________________________________________________
const char *FitCalculator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FitCalculator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FitCalculator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FitCalculator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FitCalculator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FitCalculator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FitCalculator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FitCalculator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void FitCalculator::Streamer(TBuffer &R__b)
{
   // Stream an object of class FitCalculator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FitCalculator::Class(),this);
   } else {
      R__b.WriteClassBuffer(FitCalculator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FitCalculator(void *p) {
      return  p ? new(p) ::FitCalculator : new ::FitCalculator;
   }
   static void *newArray_FitCalculator(Long_t nElements, void *p) {
      return p ? new(p) ::FitCalculator[nElements] : new ::FitCalculator[nElements];
   }
   // Wrapper around operator delete
   static void delete_FitCalculator(void *p) {
      delete ((::FitCalculator*)p);
   }
   static void deleteArray_FitCalculator(void *p) {
      delete [] ((::FitCalculator*)p);
   }
   static void destruct_FitCalculator(void *p) {
      typedef ::FitCalculator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FitCalculator

namespace {
  void TriggerDictionaryInitialization_libFitCalculator_Impl() {
    static const char* headers[] = {
"/home/nicolin/CERN/LimitSetting/src/FitCalculator.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/local/include",
"/home/nicolin/CERN/LimitSetting/bin/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libFitCalculator dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@/home/nicolin/CERN/LimitSetting/src/FitCalculator.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(implements the likelihood scan)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/nicolin/CERN/LimitSetting/src/FitCalculator.h")))  FitCalculator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libFitCalculator dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/nicolin/CERN/LimitSetting/src/FitCalculator.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"FitCalculator", payloadCode, "@",
"FitCalculator::fgIsA", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libFitCalculator",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libFitCalculator_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libFitCalculator_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libFitCalculator() {
  TriggerDictionaryInitialization_libFitCalculator_Impl();
}
