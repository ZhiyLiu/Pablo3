# This is an implementation detail for using MinimalGeodesic with the
# FindMinimalGeodesic.cmake module.  Do not include directly by name.  
# This should be included only when FindMinimalGeodesic.cmake sets 
# the MinimalGeodesic_USE_FILE variable to point here.

IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "Using MinimalGeodesic")
ENDIF(CREA_VERBOSE_CMAKE)

# Load the compiler settings used for MinimalGeodesic.
IF(MinimalGeodesic_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${MinimalGeodesic_BUILD_SETTINGS_FILE})
ENDIF(MinimalGeodesic_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use MinimalGeodesic.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MinimalGeodesic_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MinimalGeodesic_REQUIRED_CXX_FLAGS}")
SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${MinimalGeodesic_REQUIRED_LINK_FLAGS}")

# Add include directories needed to use MinimalGeodesic.
INCLUDE_DIRECTORIES(${MinimalGeodesic_INCLUDE_DIRS})

# Add link directories needed to use MinimalGeodesic.
LINK_DIRECTORIES(${MinimalGeodesic_LIBRARY_DIRS})

# Set the version 
# Already done in bbtkConfigure.h
#ADD_DEFINITIONS( -DMinimalGeodesic_VERSION="\"${MinimalGeodesic_VERSION}\"" )

# Additional use file 
IF (MinimalGeodesic_HAS_ADDITIONAL_CONFIG_FILE)
  # Include it
  INCLUDE(${MinimalGeodesic_DIR}/AdditionalUseMinimalGeodesic.cmake)
ENDIF (MinimalGeodesic_HAS_ADDITIONAL_CONFIG_FILE)
