// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__LimitCalculator

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
#include "/home/nicolin/CERN/LimitSetting/src/LimitCalculator.h"

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
   static void *new_LimitCalculator(void *p = 0);
   static void *newArray_LimitCalculator(Long_t size, void *p);
   static void delete_LimitCalculator(void *p);
   static void deleteArray_LimitCalculator(void *p);
   static void destruct_LimitCalculator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LimitCalculator*)
   {
      ::LimitCalculator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LimitCalculator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LimitCalculator", ::LimitCalculator::Class_Version(), "", 19,
                  typeid(::LimitCalculator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LimitCalculator::Dictionary, isa_proxy, 4,
                  sizeof(::LimitCalculator) );
      instance.SetNew(&new_LimitCalculator);
      instance.SetNewArray(&newArray_LimitCalculator);
      instance.SetDelete(&delete_LimitCalculator);
      instance.SetDeleteArray(&deleteArray_LimitCalculator);
      instance.SetDestructor(&destruct_LimitCalculator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LimitCalculator*)
   {
      return GenerateInitInstanceLocal((::LimitCalculator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LimitCalculator*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LimitCalculator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LimitCalculator::Class_Name()
{
   return "LimitCalculator";
}

//______________________________________________________________________________
const char *LimitCalculator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LimitCalculator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LimitCalculator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LimitCalculator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LimitCalculator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LimitCalculator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LimitCalculator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LimitCalculator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LimitCalculator::Streamer(TBuffer &R__b)
{
   // Stream an object of class LimitCalculator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LimitCalculator::Class(),this);
   } else {
      R__b.WriteClassBuffer(LimitCalculator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LimitCalculator(void *p) {
      return  p ? new(p) ::LimitCalculator : new ::LimitCalculator;
   }
   static void *newArray_LimitCalculator(Long_t nElements, void *p) {
      return p ? new(p) ::LimitCalculator[nElements] : new ::LimitCalculator[nElements];
   }
   // Wrapper around operator delete
   static void delete_LimitCalculator(void *p) {
      delete ((::LimitCalculator*)p);
   }
   static void deleteArray_LimitCalculator(void *p) {
      delete [] ((::LimitCalculator*)p);
   }
   static void destruct_LimitCalculator(void *p) {
      typedef ::LimitCalculator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LimitCalculator

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
  void TriggerDictionaryInitialization_libLimitCalculator_Impl() {
    static const char* headers[] = {
"/home/nicolin/CERN/LimitSetting/src/LimitCalculator.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/local/include",
"/home/nicolin/CERN/LimitSetting/bin/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLimitCalculator dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@/home/nicolin/CERN/LimitSetting/src/LimitCalculator.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(implements the limit calculator)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$/home/nicolin/CERN/LimitSetting/src/LimitCalculator.h")))  LimitCalculator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLimitCalculator dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/nicolin/CERN/LimitSetting/src/LimitCalculator.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"LimitCalculator", payloadCode, "@",
"LimitCalculator::fgIsA", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLimitCalculator",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLimitCalculator_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLimitCalculator_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLimitCalculator() {
  TriggerDictionaryInitialization_libLimitCalculator_Impl();
}
