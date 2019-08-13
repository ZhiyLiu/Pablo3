# - Find a library installation or build tree.
# 
# The following variables are set if MinimalGeodesic is found.  
# If MinimalGeodesic is not found, MinimalGeodesic_FOUND is set to false.
#  MinimalGeodesic_FOUND         - Set to true when MinimalGeodesic is found.
#  MinimalGeodesic_USE_FILE      - CMake file to use MinimalGeodesic.
#  MinimalGeodesic_MAJOR_VERSION - The MinimalGeodesic major version number.
#  MinimalGeodesic_MINOR_VERSION - The MinimalGeodesic minor version number 
#                       (odd non-release).
#  MinimalGeodesic_BUILD_VERSION - The MinimalGeodesic patch level 
#                       (meaningless for odd minor).
#  MinimalGeodesic_INCLUDE_DIRS  - Include directories for MinimalGeodesic
#  MinimalGeodesic_LIBRARY_DIRS  - Link directories for MinimalGeodesic libraries
#  MinimalGeodesic_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate MinimalGeodesic:
#  MinimalGeodesic_DIR  - The directory containing MinimalGeodesicConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/MinimalGeodesic directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(MinimalGeodesic_DIR_DESCRIPTION "directory containing MinimalGeodesicConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/MinimalGeodesic for an installation.")
SET(MinimalGeodesic_NOT_FOUND_MESSAGE "MinimalGeodesic not found.  Set the MinimalGeodesic_DIR cmake cache entry to the ${MinimalGeodesic_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT MinimalGeodesic_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" MinimalGeodesic_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" MinimalGeodesic_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" MinimalGeodesic_DIR_SEARCH2 "${MinimalGeodesic_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(MinimalGeodesic_DIR_SEARCH "")
  FOREACH(dir ${MinimalGeodesic_DIR_SEARCH2})
    SET(MinimalGeodesic_DIR_SEARCH ${MinimalGeodesic_DIR_SEARCH}
      ${dir}/../lib/lib64/MinimalGeodesic
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(MinimalGeodesic_DIR UseMinimalGeodesic.cmake
    # Look for an environment variable MinimalGeodesic_DIR.
    $ENV{MinimalGeodesic_DIR}

    # Look in places relative to the system executable search path.
    ${MinimalGeodesic_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/MinimalGeodesic"

    # Look in standard UNIX install locations.
    /usr/local/lib64/MinimalGeodesic
    /usr/lib64/MinimalGeodesic

    # Read from the CMakeSetup registry entries.  It is likely that
    # MinimalGeodesic will have been recently built.
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10]

    # Help the user find it if we cannot.
    DOC "The ${MinimalGeodesic_DIR_DESCRIPTION}"
  )
ENDIF(NOT MinimalGeodesic_DIR)

# If MinimalGeodesic was found, load the configuration file to get the rest of the
# settings.
IF(MinimalGeodesic_DIR)
  # Make sure the MinimalGeodesicConfig.cmake file exists in the directory provided.
  IF(EXISTS ${MinimalGeodesic_DIR}/MinimalGeodesicConfig.cmake)

    # We found MinimalGeodesic.  Load the settings.
    SET(MinimalGeodesic_FOUND 1)
    INCLUDE(${MinimalGeodesic_DIR}/MinimalGeodesicConfig.cmake)

  ENDIF(EXISTS ${MinimalGeodesic_DIR}/MinimalGeodesicConfig.cmake)
ELSE(MinimalGeodesic_DIR)
  # We did not find MinimalGeodesic.
  SET(MinimalGeodesic_FOUND 0)
ENDIF(MinimalGeodesic_DIR)

#-----------------------------------------------------------------------------
IF(NOT MinimalGeodesic_FOUND)
  # MinimalGeodesic not found, explain to the user how to specify its location.
  IF(NOT MinimalGeodesic_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${MinimalGeodesic_NOT_FOUND_MESSAGE})
  ELSE(NOT MinimalGeodesic_FIND_QUIETLY)
    IF(MinimalGeodesic_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${MinimalGeodesic_NOT_FOUND_MESSAGE})
    ENDIF(MinimalGeodesic_FIND_REQUIRED)
  ENDIF(NOT MinimalGeodesic_FIND_QUIETLY)
ENDIF(NOT MinimalGeodesic_FOUND)
