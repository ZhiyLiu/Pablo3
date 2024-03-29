#-----------------------------------------------------------------------------
#
# matchConfig.cmake - CMake configuration file for external projects.
# This file was automatically generated by the cmake macro 
# CREA_INSTALL_LIBRARY_FOR_CMAKE of the package CREA
#
# This file is configured by cmake and used by the
# Usematch.cmake module to load the lib settings 
# for an external project.

# Build tree config ?
SET(CILC_BUILD_TREE_CONFIGURATION FALSE)


IF(UNIX)
SET(GOTO_INSTALL_PREFIX /../../..)
ENDIF(UNIX)


# The match include file *RELATIVE* directories.
SET(CILC_RELATIVE_INCLUDE_PATHS "include/match")
# Compute the prefix for include and library paths
IF(CILC_BUILD_TREE_CONFIGURATION)
  # In build tree 
  # the include paths are relative to the source tree *AND* the binary tree 
  # for generated files 
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/software/sreps/Pablo2)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(match_INCLUDE_DIRS
      ${match_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
  SET(CILC_INCLUDE_PATH_PREFIX /playpen/software/sreps/Pablo2/cmake-build-debug)
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(match_INCLUDE_DIRS
      ${match_INCLUDE_DIRS}
      ${CILC_INCLUDE_PATH_PREFIX}/${path}
      )
  ENDFOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
ELSE(CILC_BUILD_TREE_CONFIGURATION)
  # In install tree 
  # the include paths are relative to install prefix 
  # On unix , GOTO_INSTALL_PREFIX allows to get back to the 
  # installation prefix from  match_DIR
  SET(CILC_INCLUDE_PATH_PREFIX ${match_DIR}${GOTO_INSTALL_PREFIX})
  # Build the *ABSOLUTE* directories
  FOREACH(path ${CILC_RELATIVE_INCLUDE_PATHS})
    SET(match_INCLUDE_DIRS
      ${match_INCLUDE_DIRS}
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
  SET(CILC_LIBRARY_PATH_PREFIX ${match_DIR}${GOTO_INSTALL_PREFIX})
ENDIF(CILC_BUILD_TREE_CONFIGURATION)
# The match library file *RELATIVE* directories.
SET(CILC_RELATIVE_LIBRARY_PATHS "lib64")
# Build the *ABSOLUTE* directories
FOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})
  SET(match_LIBRARY_DIRS
    ${match_LIBRARY_DIRS}
    ${CILC_LIBRARY_PATH_PREFIX}/${path}
    )
ENDFOREACH(path ${CILC_RELATIVE_LIBRARY_PATHS})

# Set the "prefix path"
SET(match_INSTALL_PREFIX ${CILC_LIBRARY_PATH_PREFIX})

# The C and C++ flags added by match to the cmake-configured flags.
SET(match_REQUIRED_C_FLAGS "")
SET(match_REQUIRED_CXX_FLAGS "")
SET(match_REQUIRED_LINK_FLAGS "")

# The match version 
SET(match_MAJOR_VERSION )
SET(match_MINOR_VERSION )
SET(match_BUILD_VERSION )
SET(match_VERSION ..)

# The location of the Usematch.cmake file.
SET(match_USE_FILE "${match_DIR}/Usematch.cmake")

# The build settings file.
SET(match_BUILD_SETTINGS_FILE 
  "${match_DIR}/matchBuildSettings.cmake")

# A list of all libraries for match.  Those listed here should
# automatically pull in their dependencies.
SET(match_LIBRARIES match)

# Messages
IF(CREA_VERBOSE_CMAKE)
  MESSAGE(STATUS "=======================================")
  MESSAGE(STATUS "Looking for match... found:")
  MESSAGE(STATUS "* match_DIR          = ${match_DIR}")
  MESSAGE(STATUS "* match_VERSION      = ${match_VERSION}")
  MESSAGE(STATUS "* match_USE_FILE     = ${match_USE_FILE}")

  MESSAGE(STATUS "* match_INCLUDE_DIRS = ${match_INCLUDE_DIRS}")
  MESSAGE(STATUS "* match_LIBRARY_DIRS = ${match_LIBRARY_DIRS}")
  MESSAGE(STATUS "* match_LIBRARIES    = ${match_LIBRARIES}")
ENDIF(CREA_VERBOSE_CMAKE)

# Does the library has an additional config file (user provided) ?
SET(match_HAS_ADDITIONAL_CONFIG_FILE FALSE)

IF (match_HAS_ADDITIONAL_CONFIG_FILE)
  IF(CREA_VERBOSE_CMAKE)
    MESSAGE(STATUS "Reading match additional configuration file")
  ENDIF(CREA_VERBOSE_CMAKE)
  # Include it
  INCLUDE(${match_DIR}/AdditionalmatchConfig.cmake)
ENDIF (match_HAS_ADDITIONAL_CONFIG_FILE)
