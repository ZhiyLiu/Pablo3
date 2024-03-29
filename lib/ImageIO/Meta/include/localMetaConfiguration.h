#ifdef METAIO_USE_NAMESPACE
  #undef METAIO_USE_NAMESPACE 
#endif
#ifdef METAIO_NAMESPACE
  #undef METAIO_NAMESPACE 
#endif
#ifdef METAIO_STL
  #undef METAIO_STL 
#endif
#ifdef METAIO_STREAM
  #undef METAIO_STREAM 
#endif
#ifdef METAIO_EXPORT
  #undef METAIO_EXPORT
#endif 

#include "metaIOConfig.h"

#if defined(METAIO_FOR_ITK) || !defined(METAIO_FOR_VTK)
  // ITK
  
  #define METAIO_USE_NAMESPACE  0
  #define METAIO_NAMESPACE      ITKMetaIO
  
//  #include "../itk_zlib.h"
  #include "zlib.h"
//  #include "itksys/SystemTools.hxx"
  
  #define METAIO_STL    std
  #define METAIO_STREAM std
  #define METAIO_KWSYS  itksys
  #include <iostream>
  #include <fstream>
  
  #define METAIO_EXPORT 
  
#else
// VTK

  #define METAIO_USE_NAMESPACE  1
  #define METAIO_NAMESPACE      vtkmetaio

  #include "vtkConfigure.h"
  #include "vtk_zlib.h"
  #include "vtksys/SystemTools.hxx"
  
  #ifdef VTK_NO_STD_NAMESPACE
    #define METAIO_STL
  #else
    #define METAIO_STL  std
  #endif
  
  #ifdef VTK_USE_ANSI_STDLIB
    #define METAIO_STREAM std
    #include <iostream>
    #include <fstream>
  #else
    #define METAIO_STREAM 
    #include <iostream.h>
    #include <fstream.h>
  #endif

  #define METAIO_KWSYS  vtksys
  
  #if (defined(_WIN32) || defined(WIN32)) && defined(vtkmetaio_BUILD_SHARED_LIBS)
    #ifdef vtkmetaio_EXPORTS
      #define METAIO_EXPORT __declspec(dllexport)
      #define METAIO_EXTERN
    #else
      #define METAIO_EXPORT __declspec(dllimport)
      #define METAIO_EXTERN extern
    #endif
  #else
    #define METAIO_EXPORT 
  #endif

// end VTK/ITK
#endif
