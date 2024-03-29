#-----------------------------------------------------------------------------
#
# SRepInterpolationConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# UseSRepInterpolation.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION FALSE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The SRepInterpolation include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "include/SRepInterpolation")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/software/sreps/Pablo2)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepInterpolation_INCLUDE_DIRS
      ${SRepInterpolation_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/software/sreps/Pablo2/cmake-build-debug)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepInterpolation_INCLUDE_DIRS
      ${SRepInterpolation_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  SRepInterpolation_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${SRepInterpolation_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepInterpolation_INCLUDE_DIRS
      ${SRepInterpolation_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)



# Compute the prefix for library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree
  # the library paths are relative to the binary tree 
  SET(CILC_LIBRARY_PATH_PREFIX /playpen/software/sreps/Pablo2/cmake-build-debug)
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the library paths are relative to install prefix
  SET(CILC_LIBRARY_PATH_PREFIX ${SRepInterpolation_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The SRepInterpolation library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS "lib64")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(SRepInterpolation_LIBRARY_DIRS
    ${SRepInterpolation_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(SRepInterpolation_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by SRepInterpolation to the cmake-configured flags.
SET(SRepInterpolation_REQUIRED_C_FLAGS "")
SET(SRepInterpolation_REQUIRED_CXX_FLAGS "")
SET(SRepInterpolation_REQUIRED_LINK_FLAGS "")

# The SRepInterpolation version 
SET(SRepInterpolation_MAJOR_VERSION )
SET(SRepInterpolation_MINOR_VERSION )
SET(SRepInterpolation_BUILD_VERSION )
SET(SRepInterpolation_VERSION ..)

# The location of the UseSRepInterpolation.cmake file.
SET(SRepInterpolation_USE_FILE "${SRepInterpolation_DIR}/UseSRepInterpolation.cmake")

# The build settings file.
SET(SRepInterpolation_BUILD_SETTINGS_FILE 
  "${SRepInterpolation_DIR}/SRepInterpolationBuildSettings.cmake")

# A list of all libraries for SRepInterpolation.  Those listed here should
# automatically pull in their dependencies.
SET(SRepInterpolation_LIBRARIES SRepInterpolation)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for SRepInterpolation... found:")
  MESSAGE(STATUS "* SRepInterpolation_DIR          = ${SRepInterpolation_DIR}")
  MESSAGE(STATUS "* SRepInterpolation_VERSION      = ${SRepInterpolation_VERSION}")
  MESSAGE(STATUS "* SRepInterpolation_USE_FILE     = ${SRepInterpolation_USE_FILE}")

  MESSAGE(STATUS "* SRepInterpolation_INCLUDE_DIRS = ${SRepInterpolation_INCLUDE_DIRS}")
  MESSAGE(STATUS "* SRepInterpolation_LIBRARY_DIRS = ${SRepInterpolation_LIBRARY_DIRS}")
  MESSAGE(STATUS "* SRepInterpolation_LIBRARIES    = ${SRepInterpolation_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(SRepInterpolation_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (SRepInterpolation_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading SRepInterpolation additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${SRepInterpolation_DIR}/AdditionalSRepInterpolationConfig.cmake)
ENDIF (SRepInterpolation_HAS_ADDITIONAL_CONFIG_FILE)
