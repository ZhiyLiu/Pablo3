# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /playpen/software/sreps/Pablo2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /playpen/software/sreps/Pablo2/cmake-build-debug

# Include any dependencies generated for this target.
include lib/version/CMakeFiles/version.dir/depend.make

# Include the progress variables for this target.
include lib/version/CMakeFiles/version.dir/progress.make

# Include the compile flags for this target's objects.
include lib/version/CMakeFiles/version.dir/flags.make

lib/version/CMakeFiles/version.dir/src/version.cpp.o: lib/version/CMakeFiles/version.dir/flags.make
lib/version/CMakeFiles/version.dir/src/version.cpp.o: ../lib/version/src/version.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/software/sreps/Pablo2/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/version/CMakeFiles/version.dir/src/version.cpp.o"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/version.dir/src/version.cpp.o -c /playpen/software/sreps/Pablo2/lib/version/src/version.cpp

lib/version/CMakeFiles/version.dir/src/version.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/version.dir/src/version.cpp.i"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/software/sreps/Pablo2/lib/version/src/version.cpp > CMakeFiles/version.dir/src/version.cpp.i

lib/version/CMakeFiles/version.dir/src/version.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/version.dir/src/version.cpp.s"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/software/sreps/Pablo2/lib/version/src/version.cpp -o CMakeFiles/version.dir/src/version.cpp.s

lib/version/CMakeFiles/version.dir/src/version.cpp.o.requires:
.PHONY : lib/version/CMakeFiles/version.dir/src/version.cpp.o.requires

lib/version/CMakeFiles/version.dir/src/version.cpp.o.provides: lib/version/CMakeFiles/version.dir/src/version.cpp.o.requires
	$(MAKE) -f lib/version/CMakeFiles/version.dir/build.make lib/version/CMakeFiles/version.dir/src/version.cpp.o.provides.build
.PHONY : lib/version/CMakeFiles/version.dir/src/version.cpp.o.provides

lib/version/CMakeFiles/version.dir/src/version.cpp.o.provides.build: lib/version/CMakeFiles/version.dir/src/version.cpp.o

# Object files for target version
version_OBJECTS = \
"CMakeFiles/version.dir/src/version.cpp.o"

# External object files for target version
version_EXTERNAL_OBJECTS =

libversion.a: lib/version/CMakeFiles/version.dir/src/version.cpp.o
libversion.a: lib/version/CMakeFiles/version.dir/build.make
libversion.a: lib/version/CMakeFiles/version.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../libversion.a"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version && $(CMAKE_COMMAND) -P CMakeFiles/version.dir/cmake_clean_target.cmake
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/version.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/version/CMakeFiles/version.dir/build: libversion.a
.PHONY : lib/version/CMakeFiles/version.dir/build

lib/version/CMakeFiles/version.dir/requires: lib/version/CMakeFiles/version.dir/src/version.cpp.o.requires
.PHONY : lib/version/CMakeFiles/version.dir/requires

lib/version/CMakeFiles/version.dir/clean:
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version && $(CMAKE_COMMAND) -P CMakeFiles/version.dir/cmake_clean.cmake
.PHONY : lib/version/CMakeFiles/version.dir/clean

lib/version/CMakeFiles/version.dir/depend:
	cd /playpen/software/sreps/Pablo2/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /playpen/software/sreps/Pablo2 /playpen/software/sreps/Pablo2/lib/version /playpen/software/sreps/Pablo2/cmake-build-debug /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version /playpen/software/sreps/Pablo2/cmake-build-debug/lib/version/CMakeFiles/version.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/version/CMakeFiles/version.dir/depend

