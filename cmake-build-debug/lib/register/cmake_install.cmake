# Install script for directory: /playpen/software/sreps/Pablo2/lib/register

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/register" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/register/InstallLibraryForCMake_tmp/Useregister.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/register" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/register/InstallLibraryForCMake_tmp/registerConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/register" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/register/InstallLibraryForCMake_tmp/registerBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findregister.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/register/InstallLibraryForCMake_tmp/Findregister.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/register" TYPE FILE FILES
    "/playpen/software/sreps/Pablo2/lib/register/include/FunctionExplorer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSRepProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/LandmarkDeformation.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DDeformationPGAProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSRepOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DRegistrationProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSubfigureProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DFigureElongater.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSpokeProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DVoxelOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/Trackball.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DtoPovray.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/PCA.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSubfigureTransformation.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DMainFigResiduePGAProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/Anastruct.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSubfigureOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/P3DControl.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/OptimizerBase.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DRegistrationOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/P3DControlBkp1.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DRegistrationPGAOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DMainFigureProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSimilarityPGAProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DMainFigResidueCPNSProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DBoundingSphere.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/OctTreeTauBand.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSimilarityCPNSOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/ImagePlanes.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSimilarityElongationProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSimilarityPGAOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DDeformationProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DRegistrationSimilarityProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DPrimitiveCorrector.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/FitUnlabeledPoints.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSimilarityElongationOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DTubeSimilarityPGAProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DAdaptiveRegistrationPGAOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/P3DUndoList.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DRegistrationSimilarityOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DDeformationOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/DistanceToPointSetFunction.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DMainFigureOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DAdaptiveRegistrationPGAProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSpokeOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DSimilarityCPNSProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DMainFigResidueCPNSOptimizer.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/MakeOptimizationVideo.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DRegistrationPGAProblem.h"
    "/playpen/software/sreps/Pablo2/lib/register/include/M3DMainFigResiduePGAOptimizer.h"
    "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/register/register_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/libregister.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
