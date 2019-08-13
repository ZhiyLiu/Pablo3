#ifndef __MinimalGeodesic_EXPORT_h_INCLUDED__
#define __MinimalGeodesic_EXPORT_h_INCLUDED__

// Automatically generated file which defines 
// the symbols MinimalGeodesic_EXPORT and MinimalGeodesic_CDECL
// to be used for the definition of classes or functions
// which must be exported when the lib is built as a shared lib on Windows
// and imported when the shared lib is used by another program

#if defined(_WIN32) && defined (MinimalGeodesic_BUILD_SHARED)
  #ifdef MinimalGeodesic_EXPORT_SYMBOLS
    #define MinimalGeodesic_EXPORT __declspec( dllexport )
  #else
    #define MinimalGeodesic_EXPORT __declspec( dllimport )
  #endif
  #define MinimalGeodesic_CDECL __cdecl
#else
  #define MinimalGeodesic_EXPORT
  #define MinimalGeodesic_CDECL
#endif // defined(_WIN32)

#endif

