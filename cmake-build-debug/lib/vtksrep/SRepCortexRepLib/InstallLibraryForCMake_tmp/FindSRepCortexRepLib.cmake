# - Find a library installation or build tree.
# 
# The following variables are set if SRepCortexRepLib is found.  
# If SRepCortexRepLib is not found, SRepCortexRepLib_FOUND is set to false.
#  SRepCortexRepLib_FOUND         - Set to true when SRepCortexRepLib is found.
#  SRepCortexRepLib_USE_FILE      - CMake file to use SRepCortexRepLib.
#  SRepCortexRepLib_MAJOR_VERSION - The SRepCortexRepLib major version number.
#  SRepCortexRepLib_MINOR_VERSION - The SRepCortexRepLib minor version number 
#                       (odd non-release).
#  SRepCortexRepLib_BUILD_VERSION - The SRepCortexRepLib patch level 
#                       (meaningless for odd minor).
#  SRepCortexRepLib_INCLUDE_DIRS  - Include directories for SRepCortexRepLib
#  SRepCortexRepLib_LIBRARY_DIRS  - Link directories for SRepCortexRepLib libraries
#  SRepCortexRepLib_LIBRARIES     - List of libraries to link against
#
# The following cache entries must be set by the user to locate SRepCortexRepLib:
#  SRepCortexRepLib_DIR  - The directory containing SRepCortexRepLibConfig.cmake.  
#             This is either the root of the build tree,
#             or the lib/SRepCortexRepLib directory.  This is the 
#             only cache entry.


# Construct consitent error messages for use below.
SET(SRepCortexRepLib_DIR_DESCRIPTION "directory containing SRepCortexRepLibConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/SRepCortexRepLib for an installation.")
SET(SRepCortexRepLib_NOT_FOUND_MESSAGE "SRepCortexRepLib not found.  Set the SRepCortexRepLib_DIR cmake cache entry to the ${SRepCortexRepLib_DIR_DESCRIPTION}")

# Search only if the location is not already known.
IF(NOT SRepCortexRepLib_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" SRepCortexRepLib_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" SRepCortexRepLib_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" SRepCortexRepLib_DIR_SEARCH2 "${SRepCortexRepLib_DIR_SEARCH1}")

  # Construct a set of paths relative to the system search path.
  SET(SRepCortexRepLib_DIR_SEARCH "")
  FOREACH(dir ${SRepCortexRepLib_DIR_SEARCH2})
    SET(SRepCortexRepLib_DIR_SEARCH ${SRepCortexRepLib_DIR_SEARCH}
      ${dir}/../lib/lib64/SRepCortexRepLib
      )
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(SRepCortexRepLib_DIR UseSRepCortexRepLib.cmake
    # Look for an environment variable SRepCortexRepLib_DIR.
    $ENV{SRepCortexRepLib_DIR}

    # Look in places relative to the system executable search path.
    ${SRepCortexRepLib_DIR_SEARCH}

    # Look in standard WIN install locations.
    "$ENV{ProgramFiles}/SRepCortexRepLib"

    # Look in standard UNIX install locations.
    /usr/local/lib64/SRepCortexRepLib
    /usr/lib64/SRepCortexRepLib

    # Read from the CMakeSetup registry entries.  It is likely that
    # SRepCortexRepLib will have been recently built.
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
    DOC "The ${SRepCortexRepLib_DIR_DESCRIPTION}"
  )
ENDIF(NOT SRepCortexRepLib_DIR)

# If SRepCortexRepLib was found, load the configuration file to get the rest of the
# settings.
IF(SRepCortexRepLib_DIR)
  # Make sure the SRepCortexRepLibConfig.cmake file exists in the directory provided.
  IF(EXISTS ${SRepCortexRepLib_DIR}/SRepCortexRepLibConfig.cmake)

    # We found SRepCortexRepLib.  Load the settings.
    SET(SRepCortexRepLib_FOUND 1)
    INCLUDE(${SRepCortexRepLib_DIR}/SRepCortexRepLibConfig.cmake)

  ENDIF(EXISTS ${SRepCortexRepLib_DIR}/SRepCortexRepLibConfig.cmake)
ELSE(SRepCortexRepLib_DIR)
  # We did not find SRepCortexRepLib.
  SET(SRepCortexRepLib_FOUND 0)
ENDIF(SRepCortexRepLib_DIR)

#-----------------------------------------------------------------------------
IF(NOT SRepCortexRepLib_FOUND)
  # SRepCortexRepLib not found, explain to the user how to specify its location.
  IF(NOT SRepCortexRepLib_FIND_QUIETLY)
    MESSAGE(FATAL_ERROR ${SRepCortexRepLib_NOT_FOUND_MESSAGE})
  ELSE(NOT SRepCortexRepLib_FIND_QUIETLY)
    IF(SRepCortexRepLib_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR ${SRepCortexRepLib_NOT_FOUND_MESSAGE})
    ENDIF(SRepCortexRepLib_FIND_REQUIRED)
  ENDIF(NOT SRepCortexRepLib_FIND_QUIETLY)
ENDIF(NOT SRepCortexRepLib_FOUND)
