# Install script for directory: /playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/SRepCortexRepLib" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/vtksrep/SRepCortexRepLib/InstallLibraryForCMake_tmp/UseSRepCortexRepLib.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/SRepCortexRepLib" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/vtksrep/SRepCortexRepLib/InstallLibraryForCMake_tmp/SRepCortexRepLibConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/SRepCortexRepLib" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/vtksrep/SRepCortexRepLib/InstallLibraryForCMake_tmp/SRepCortexRepLibBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/FindSRepCortexRepLib.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/vtksrep/SRepCortexRepLib/InstallLibraryForCMake_tmp/FindSRepCortexRepLib.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SRepCortexRepLib" TYPE FILE FILES
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/optimizeCortexSpokes.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/optimizeCortexAtoms.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/optimizeSRepSampling.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/vtkoptimizeCortexAtomsDis.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/vtkoptimizeCortexAtomsCurv.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/optimizeCortexAtomsDis.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/vtkoptimizeSRepSampling.h"
    "/playpen/software/sreps/Pablo2/lib/vtksrep/SRepCortexRepLib/vtkoptimizeCortexSpokes.h"
    "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/vtksrep/SRepCortexRepLib/SRepCortexRepLib_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/libSRepCortexRepLib.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

