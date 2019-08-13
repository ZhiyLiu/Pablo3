#ifndef __SRepCortexRepLib_EXPORT_h_INCLUDED__
#define __SRepCortexRepLib_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols SRepCortexRepLib_EXPORT and SRepCortexRepLib_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (SRepCortexRepLib_BUILD_SHARED)
  #ifdef SRepCortexRepLib_EXPORT_SYMBOLS
    #define SRepCortexRepLib_EXPORT __declspec( dllexport )
  #else
    #define SRepCortexRepLib_EXPORT __declspec( dllimport )
  #endif
  #define SRepCortexRepLib_CDECL __cdecl
#else
  #define SRepCortexRepLib_EXPORT
  #define SRepCortexRepLib_CDECL
#endif // defined(_WIN32)

#endif

