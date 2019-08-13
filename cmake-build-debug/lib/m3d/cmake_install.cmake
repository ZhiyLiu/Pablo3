# Install script for directory: /playpen/software/sreps/Pablo2/lib/m3d

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/m3d" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/m3d/InstallLibraryForCMake_tmp/Usem3d.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/m3d" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/m3d/InstallLibraryForCMake_tmp/m3dConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/m3d" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/m3d/InstallLibraryForCMake_tmp/m3dBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findm3d.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/m3d/InstallLibraryForCMake_tmp/Findm3d.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/m3d" TYPE FILE FILES
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DEndPrimitive.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/PluncMatrixFile.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Hermite.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DFigure.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DFigureStats.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/SubdivBoundary.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Mathdefs.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DQuadFigure.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/RigidShapeSpace.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DQuadFigureRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Quat.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/ScaleTransShapeSpace.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/SimilarityTransform3D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/RotationAndScaleShapeSpace.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DFigureRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DTubeEndPrimitive.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DTubePrimitive.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DObjectFile.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DAtomPredictorTube.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/pablo_version.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/TileSet.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DPNSTransform.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DObjectRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DPrimitiveRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/NormalAlignedShapeSpace.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/BYU.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/WorldSystem.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/DQFImage.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/utility.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DCPNSStats.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DFigurePredictor.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DSpoke_bkp_0.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DPGAStats.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DPGA.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/InterfiguralConstraints.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DPrimitive.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/ControlParmsAccess.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Bezier1D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DQuadEndPrimitive.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/RAWImageFile.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/GeodesicDistanceFunction.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Image3D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/GaussianBlur3D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Bezier2D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DQuadInterpolator.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DPGAPrimitiveStats.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DTubeFigure.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/ImageResample3D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DAtomPredictorQuad.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/control_parms_defaults.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DAtomPredictor.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Vector3D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DTubeFigureRenderer.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DInterpolater.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DQuadPrimitive.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Matrix4D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Vector2D.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/Geodesic.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DFigureTreeNode.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DObject.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/ControlParms.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/GeodesicSym.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DDisplayGlobals.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/ControlFile.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/IntervalTimer.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/VectorND.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DSpoke.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DQuadInterpolater.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/SimTransShapeSpace.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/snapshot.h"
    "/playpen/software/sreps/Pablo2/lib/m3d/include/M3DSpoke_bkp_1.h"
    "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/m3d/m3d_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/libm3d.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

