#ifndef __SRepInterpolation_EXPORT_h_INCLUDED__
#define __SRepInterpolation_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols SRepInterpolation_EXPORT and SRepInterpolation_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (SRepInterpolation_BUILD_SHARED)
  #ifdef SRepInterpolation_EXPORT_SYMBOLS
    #define SRepInterpolation_EXPORT __declspec( dllexport )
  #else
    #define SRepInterpolation_EXPORT __declspec( dllimport )
  #endif
  #define SRepInterpolation_CDECL __cdecl
#else
  #define SRepInterpolation_EXPORT
  #define SRepInterpolation_CDECL
#endif // defined(_WIN32)

#endif

