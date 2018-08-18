// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__PValueCalculator

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
#include "/home/nicolin/CERN/LimitSetting/src/PValueCalculator.h"

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
   static void *new_PValueCalculator(void *p = 0);
   static void *newArray_PValueCalculator(Long_t size, void *p);
   static void delete_PValueCalculator(void *p);
   static void deleteArray_PValueCalculator(void *p);
   static void destruct_PValueCalculator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PValueCalculator*)
   {
      ::PValueCalculator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PValueCalculator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PValueCalculator", ::PValueCalculator::Class_Version(), "", 23,
                  typeid(::PValueCalculator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PValueCalculator::Dictionary, isa_proxy, 4,
                  sizeof(::PValueCalculator) );
      instance.SetNew(&new_PValueCalculator);
      instance.SetNewArray(&newArray_PValueCalculator);
      instance.SetDelete(&delete_PValueCalculator);
      instance.SetDeleteArray(&deleteArray_PValueCalculator);
      instance.SetDestructor(&destruct_PValueCalculator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PValueCalculator*)
   {
      return GenerateInitInstanceLocal((::PValueCalculator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PValueCalculator*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PValueCalculator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PValueCalculator::Class_Name()
{
   return "PValueCalculator";
}

//______________________________________________________________________________
const char *PValueCalculator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PValueCalculator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PValueCalculator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PValueCalculator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PValueCalculator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PValueCalculator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PValueCalculator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PValueCalculator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PValueCalculator::Streamer(TBuffer &R__b)
{
   // Stream an object of class PValueCalculator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PValueCalculator::Class(),this);
   } else {
      R__b.WriteClassBuffer(PValueCalculator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PValueCalculator(void *p) {
      return  p ? new(p) ::PValueCalculator : new ::PValueCalculator;
   }
   static void *newArray_PValueCalculator(Long_t nElements, void *p) {
      return p ? new(p) ::PValueCalculator[nElements] : new ::PValueCalculator[nElements];
   }
   // Wrapper around operator delete
   static void delete_PValueCalculator(void *p) {
      delete ((::PValueCalculator*)p);
   }
   static void deleteArray_PValueCalculator(void *p) {
      delete [] ((::PValueCalculator*)p);
   }
   static void destruct_PValueCalculator(void *p) {
      typedef ::PValueCalculator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PValueCalculator

namespace {
  void TriggerDictionaryInitialization_libPValueCalculator_Impl() {
    static const char* headers[] = {
"/home/nicolin/CERN/LimitSetting/src/PValueCalculator.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/local/include",
"/home/nicolin/CERN/LimitSetting/bin/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libPValueCalculator dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@/home/nicolin/CERN/LimitSetting/src/PValueCalculator.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(implements the p-value calculator)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/nicolin/CERN/LimitSetting/src/PValueCalculator.h")))  PValueCalculator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libPValueCalculator dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/nicolin/CERN/LimitSetting/src/PValueCalculator.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"PValueCalculator", payloadCode, "@",
"PValueCalculator::fgIsA", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libPValueCalculator",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libPValueCalculator_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libPValueCalculator_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libPValueCalculator() {
  TriggerDictionaryInitialization_libPValueCalculator_Impl();
}
