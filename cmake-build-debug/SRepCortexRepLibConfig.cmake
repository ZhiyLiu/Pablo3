#-----------------------------------------------------------------------------
#
# SRepCortexRepLibConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# UseSRepCortexRepLib.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION TRUE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The SRepCortexRepLib include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "lib/vtksrep/SRepCortexRepLib")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/software/sreps/Pablo2)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepCortexRepLib_INCLUDE_DIRS
      ${SRepCortexRepLib_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/software/sreps/Pablo2/cmake-build-debug)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepCortexRepLib_INCLUDE_DIRS
      ${SRepCortexRepLib_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  SRepCortexRepLib_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${SRepCortexRepLib_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(SRepCortexRepLib_INCLUDE_DIRS
      ${SRepCortexRepLib_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${SRepCortexRepLib_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The SRepCortexRepLib library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS ".")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(SRepCortexRepLib_LIBRARY_DIRS
    ${SRepCortexRepLib_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(SRepCortexRepLib_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by SRepCortexRepLib to the cmake-configured flags.
SET(SRepCortexRepLib_REQUIRED_C_FLAGS "")
SET(SRepCortexRepLib_REQUIRED_CXX_FLAGS "")
SET(SRepCortexRepLib_REQUIRED_LINK_FLAGS "")

# The SRepCortexRepLib version 
SET(SRepCortexRepLib_MAJOR_VERSION )
SET(SRepCortexRepLib_MINOR_VERSION )
SET(SRepCortexRepLib_BUILD_VERSION )
SET(SRepCortexRepLib_VERSION ..)

# The location of the UseSRepCortexRepLib.cmake file.
SET(SRepCortexRepLib_USE_FILE "${SRepCortexRepLib_DIR}/UseSRepCortexRepLib.cmake")

# The build settings file.
SET(SRepCortexRepLib_BUILD_SETTINGS_FILE 
  "${SRepCortexRepLib_DIR}/SRepCortexRepLibBuildSettings.cmake")

# A list of all libraries for SRepCortexRepLib.  Those listed here should
# automatically pull in their dependencies.
SET(SRepCortexRepLib_LIBRARIES SRepCortexRepLib)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for SRepCortexRepLib... found:")
  MESSAGE(STATUS "* SRepCortexRepLib_DIR          = ${SRepCortexRepLib_DIR}")
  MESSAGE(STATUS "* SRepCortexRepLib_VERSION      = ${SRepCortexRepLib_VERSION}")
  MESSAGE(STATUS "* SRepCortexRepLib_USE_FILE     = ${SRepCortexRepLib_USE_FILE}")

  MESSAGE(STATUS "* SRepCortexRepLib_INCLUDE_DIRS = ${SRepCortexRepLib_INCLUDE_DIRS}")
  MESSAGE(STATUS "* SRepCortexRepLib_LIBRARY_DIRS = ${SRepCortexRepLib_LIBRARY_DIRS}")
  MESSAGE(STATUS "* SRepCortexRepLib_LIBRARIES    = ${SRepCortexRepLib_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(SRepCortexRepLib_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (SRepCortexRepLib_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading SRepCortexRepLib additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${SRepCortexRepLib_DIR}/AdditionalSRepCortexRepLibConfig.cmake)
ENDIF (SRepCortexRepLib_HAS_ADDITIONAL_CONFIG_FILE)
