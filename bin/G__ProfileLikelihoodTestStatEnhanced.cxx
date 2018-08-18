// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__ProfileLikelihoodTestStatEnhanced

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
#include "/home/nicolin/CERN/LimitSetting/src/ProfileLikelihoodTestStatEnhanced.h"

// Header files passed via #pragma extra_include

namespace RooStats {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *RooStats_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("RooStats", 0 /*version*/, "RooStats/TestStatistic.h", 22,
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
   static void *new_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p = 0);
   static void *newArray_RooStatscLcLProfileLikelihoodTestStatEnhanced(Long_t size, void *p);
   static void delete_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p);
   static void deleteArray_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p);
   static void destruct_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooStats::ProfileLikelihoodTestStatEnhanced*)
   {
      ::RooStats::ProfileLikelihoodTestStatEnhanced *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooStats::ProfileLikelihoodTestStatEnhanced >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooStats::ProfileLikelihoodTestStatEnhanced", ::RooStats::ProfileLikelihoodTestStatEnhanced::Class_Version(), "", 52,
                  typeid(::RooStats::ProfileLikelihoodTestStatEnhanced), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooStats::ProfileLikelihoodTestStatEnhanced::Dictionary, isa_proxy, 4,
                  sizeof(::RooStats::ProfileLikelihoodTestStatEnhanced) );
      instance.SetNew(&new_RooStatscLcLProfileLikelihoodTestStatEnhanced);
      instance.SetNewArray(&newArray_RooStatscLcLProfileLikelihoodTestStatEnhanced);
      instance.SetDelete(&delete_RooStatscLcLProfileLikelihoodTestStatEnhanced);
      instance.SetDeleteArray(&deleteArray_RooStatscLcLProfileLikelihoodTestStatEnhanced);
      instance.SetDestructor(&destruct_RooStatscLcLProfileLikelihoodTestStatEnhanced);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooStats::ProfileLikelihoodTestStatEnhanced*)
   {
      return GenerateInitInstanceLocal((::RooStats::ProfileLikelihoodTestStatEnhanced*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooStats::ProfileLikelihoodTestStatEnhanced*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace RooStats {
//______________________________________________________________________________
atomic_TClass_ptr ProfileLikelihoodTestStatEnhanced::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ProfileLikelihoodTestStatEnhanced::Class_Name()
{
   return "RooStats::ProfileLikelihoodTestStatEnhanced";
}

//______________________________________________________________________________
const char *ProfileLikelihoodTestStatEnhanced::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooStats::ProfileLikelihoodTestStatEnhanced*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ProfileLikelihoodTestStatEnhanced::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooStats::ProfileLikelihoodTestStatEnhanced*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ProfileLikelihoodTestStatEnhanced::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooStats::ProfileLikelihoodTestStatEnhanced*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ProfileLikelihoodTestStatEnhanced::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooStats::ProfileLikelihoodTestStatEnhanced*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RooStats
namespace RooStats {
//______________________________________________________________________________
void ProfileLikelihoodTestStatEnhanced::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooStats::ProfileLikelihoodTestStatEnhanced.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RooStats::ProfileLikelihoodTestStatEnhanced::Class(),this);
   } else {
      R__b.WriteClassBuffer(RooStats::ProfileLikelihoodTestStatEnhanced::Class(),this);
   }
}

} // namespace RooStats
namespace ROOT {
   // Wrappers around operator new
   static void *new_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RooStats::ProfileLikelihoodTestStatEnhanced : new ::RooStats::ProfileLikelihoodTestStatEnhanced;
   }
   static void *newArray_RooStatscLcLProfileLikelihoodTestStatEnhanced(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RooStats::ProfileLikelihoodTestStatEnhanced[nElements] : new ::RooStats::ProfileLikelihoodTestStatEnhanced[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p) {
      delete ((::RooStats::ProfileLikelihoodTestStatEnhanced*)p);
   }
   static void deleteArray_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p) {
      delete [] ((::RooStats::ProfileLikelihoodTestStatEnhanced*)p);
   }
   static void destruct_RooStatscLcLProfileLikelihoodTestStatEnhanced(void *p) {
      typedef ::RooStats::ProfileLikelihoodTestStatEnhanced current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RooStats::ProfileLikelihoodTestStatEnhanced

namespace {
  void TriggerDictionaryInitialization_libProfileLikelihoodTestStatEnhanced_Impl() {
    static const char* headers[] = {
"/home/nicolin/CERN/LimitSetting/src/ProfileLikelihoodTestStatEnhanced.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/local/include",
"/home/nicolin/CERN/LimitSetting/bin/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libProfileLikelihoodTestStatEnhanced dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace RooStats{class __attribute__((annotate(R"ATTRDUMP(file_name@@@/home/nicolin/CERN/LimitSetting/src/ProfileLikelihoodTestStatEnhanced.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(implements the enhanced profile likelihood ratio as a test statistic to be used with several tools)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/nicolin/CERN/LimitSetting/src/ProfileLikelihoodTestStatEnhanced.h")))  ProfileLikelihoodTestStatEnhanced;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libProfileLikelihoodTestStatEnhanced dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/nicolin/CERN/LimitSetting/src/ProfileLikelihoodTestStatEnhanced.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"RooStats::ProfileLikelihoodTestStatEnhanced", payloadCode, "@",
"RooStats::ProfileLikelihoodTestStatEnhanced::LimitType", payloadCode, "@",
"RooStats::ProfileLikelihoodTestStatEnhanced::OptimizeBits", payloadCode, "@",
"RooStats::ProfileLikelihoodTestStatEnhanced::default_optimize", payloadCode, "@",
"RooStats::ProfileLikelihoodTestStatEnhanced::fgAlwaysReuseNll", payloadCode, "@",
"RooStats::ProfileLikelihoodTestStatEnhanced::fgIsA", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libProfileLikelihoodTestStatEnhanced",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libProfileLikelihoodTestStatEnhanced_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libProfileLikelihoodTestStatEnhanced_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libProfileLikelihoodTestStatEnhanced() {
  TriggerDictionaryInitialization_libProfileLikelihoodTestStatEnhanced_Impl();
}
