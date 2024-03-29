# Install script for directory: /playpen/software/sreps/Pablo2/lib/seurat

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Debug")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/seurat" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/seurat/InstallLibraryForCMake_tmp/Useseurat.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/seurat" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/seurat/InstallLibraryForCMake_tmp/seuratConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/seurat" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/seurat/InstallLibraryForCMake_tmp/seuratBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findseurat.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/seurat/InstallLibraryForCMake_tmp/Findseurat.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seurat" TYPE FILE FILES
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Shapedepend.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/MyQueue.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/HanInterpolation.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/M3DBlendedRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/PseudoSet.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Pointlist_serverB.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/M3DObjectSurfaceVisualizer.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Conjgrad2.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/CCSubdivsurf.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Zerofinder.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/SelectedPartialFigures.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/TileSetRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Diatomgrid.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/renderDefinitions.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Intersection.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Diatom.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/LinAlg.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Xferlist.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/M3DObjectSurfaceRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/SurfaceColorMap.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Plist_subdivcomp.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Tritri.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Samplestats.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Mesh.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/QuadMesh.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/MyList.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Subdivsurf.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/TubeMesh.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Shapeheader.h"
    "/playpen/software/sreps/Pablo2/lib/seurat/include/Pointlist_server2.h"
    "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/seurat/seurat_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/libseurat.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

