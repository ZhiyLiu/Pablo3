# Install script for directory: /playpen/software/sreps/Pablo2/lib/match

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/match" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/match/InstallLibraryForCMake_tmp/Usematch.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/match" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/match/InstallLibraryForCMake_tmp/matchConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/match" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/match/InstallLibraryForCMake_tmp/matchBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/Findmatch.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/match/InstallLibraryForCMake_tmp/Findmatch.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/match" TYPE FILE FILES
    "/playpen/software/sreps/Pablo2/lib/match/include/gpTuning.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/DQFMatch.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/deltas.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/MethodOfMoments.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/TemplateProfiles.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/MaskFile.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/SimpleMaskFile.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/DistanceMap3D.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/bpTuning.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/Tuning.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/DistanceVectorList.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/ImageDistanceMap.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/SurfacePatchEnsemble.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/Match.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/Danielsson.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/ObjectRelativeSampling.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/MatchUtility.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/ccl.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/Curvature.h"
    "/playpen/software/sreps/Pablo2/lib/match/include/Mask.h"
    "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/match/match_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/libmatch.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

