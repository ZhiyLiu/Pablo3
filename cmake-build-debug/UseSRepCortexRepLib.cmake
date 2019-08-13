# This is an implementation detail for using SRepCortexRepLib with the
# FindSRepCortexRepLib.cmake module.  Do not include directly by name.  
# This should be included only when FindSRepCortexRepLib.cmake sets 
# the SRepCortexRepLib_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using SRepCortexRepLib")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for SRepCortexRepLib.
IF(SRepCortexRepLib_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${SRepCortexRepLib_BUILD_SETTINGS_FILE})
ENDIF(SRepCortexRepLib_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use SRepCortexRepLib.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SRepCortexRepLib_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SRepCortexRepLib_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${SRepCortexRepLib_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use SRepCortexRepLib.
INCLUDE_DIRECTORIES(${SRepCortexRepLib_INCLUDE_DIRS})

# Add link directories needed to use SRepCortexRepLib.
LINK_DIRECTORIES(${SRepCortexRepLib_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DSRepCortexRepLib_VERSION="\"${SRepCortexRepLib_VERSION}\"" )

# Additional use file 
IF (SRepCortexRepLib_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${SRepCortexRepLib_DIR}/AdditionalUseSRepCortexRepLib.cmake)
ENDIF (SRepCortexRepLib_HAS_ADDITIONAL_CONFIG_FILE)
