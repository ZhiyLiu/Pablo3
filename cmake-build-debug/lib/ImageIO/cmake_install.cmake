# Install script for directory: /playpen/software/sreps/Pablo2/lib/ImageIO

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
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/ImageIO" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/ImageIO/InstallLibraryForCMake_tmp/UseImageIO.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/ImageIO" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/ImageIO/InstallLibraryForCMake_tmp/ImageIOConfig.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/ImageIO" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/ImageIO/InstallLibraryForCMake_tmp/ImageIOBuildSettings.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/cmake/Modules/FindImageIO.cmake")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/share/cmake/Modules" TYPE FILE FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/ImageIO/InstallLibraryForCMake_tmp/FindImageIO.cmake")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ImageIO" TYPE FILE FILES
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Gipl/include/ipimage.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Gipl/include/giplio.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Gipl/include/ipmatrix.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Gipl/include/misc.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Gipl/include/gipl_header.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/include/const.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/include/ImageStruct.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/include/AllImageIO.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/include/macros.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/include/ImageIO.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaIOConfig.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaUtils.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaImageUtils.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaEvent.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaImageTypes.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/localMetaConfiguration.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaObject.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaImage.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Meta/include/metaTypes.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Analyze/include/Analyze.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Analyze/include/mayo_analyze.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/Analyze/include/analyze_io.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/utilities/include/BinaryIO.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/utilities/include/BasicException.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/utilities/include/Debugging.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/getopt.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/libplan.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/xdr_ll_planio.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/xdr.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/extbm.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/gen.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/unistd.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/types.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/plan_sys.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/plan_file_io.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/plan_strings.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/libmisc.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/plan_xdr.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/libplanio.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/plan_config.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/strings.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/plan_im.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/rpc_winnt/xdr.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/rpc_winnt/types.h"
    "/playpen/software/sreps/Pablo2/lib/ImageIO/PlanIm/include/libbrachy.h"
    "/playpen/software/sreps/Pablo2/cmake-build-debug/lib/ImageIO/ImageIO_EXPORT.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/playpen/software/sreps/Pablo2/cmake-build-debug/libImageIO.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

