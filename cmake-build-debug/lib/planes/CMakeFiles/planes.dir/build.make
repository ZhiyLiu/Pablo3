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
include lib/planes/CMakeFiles/planes.dir/depend.make

# Include the progress variables for this target.
include lib/planes/CMakeFiles/planes.dir/progress.make

# Include the compile flags for this target's objects.
include lib/planes/CMakeFiles/planes.dir/flags.make

lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o: lib/planes/CMakeFiles/planes.dir/flags.make
lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o: ../lib/planes/src/CutPlanes.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/software/sreps/Pablo2/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/planes.dir/src/CutPlanes.cpp.o -c /playpen/software/sreps/Pablo2/lib/planes/src/CutPlanes.cpp

lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/planes.dir/src/CutPlanes.cpp.i"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/software/sreps/Pablo2/lib/planes/src/CutPlanes.cpp > CMakeFiles/planes.dir/src/CutPlanes.cpp.i

lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/planes.dir/src/CutPlanes.cpp.s"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/software/sreps/Pablo2/lib/planes/src/CutPlanes.cpp -o CMakeFiles/planes.dir/src/CutPlanes.cpp.s

lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.requires:
.PHONY : lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.requires

lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.provides: lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.requires
	$(MAKE) -f lib/planes/CMakeFiles/planes.dir/build.make lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.provides.build
.PHONY : lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.provides

lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.provides.build: lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o

lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o: lib/planes/CMakeFiles/planes.dir/flags.make
lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o: ../lib/planes/src/drawAtomCut.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/software/sreps/Pablo2/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/planes.dir/src/drawAtomCut.cpp.o -c /playpen/software/sreps/Pablo2/lib/planes/src/drawAtomCut.cpp

lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/planes.dir/src/drawAtomCut.cpp.i"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/software/sreps/Pablo2/lib/planes/src/drawAtomCut.cpp > CMakeFiles/planes.dir/src/drawAtomCut.cpp.i

lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/planes.dir/src/drawAtomCut.cpp.s"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/software/sreps/Pablo2/lib/planes/src/drawAtomCut.cpp -o CMakeFiles/planes.dir/src/drawAtomCut.cpp.s

lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.requires:
.PHONY : lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.requires

lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.provides: lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.requires
	$(MAKE) -f lib/planes/CMakeFiles/planes.dir/build.make lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.provides.build
.PHONY : lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.provides

lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.provides.build: lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o

lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o: lib/planes/CMakeFiles/planes.dir/flags.make
lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o: ../lib/planes/src/P3DCutPlaneView.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/software/sreps/Pablo2/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o -c /playpen/software/sreps/Pablo2/lib/planes/src/P3DCutPlaneView.cpp

lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.i"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/software/sreps/Pablo2/lib/planes/src/P3DCutPlaneView.cpp > CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.i

lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.s"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/software/sreps/Pablo2/lib/planes/src/P3DCutPlaneView.cpp -o CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.s

lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.requires:
.PHONY : lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.requires

lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.provides: lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.requires
	$(MAKE) -f lib/planes/CMakeFiles/planes.dir/build.make lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.provides.build
.PHONY : lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.provides

lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.provides.build: lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o

lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o: lib/planes/CMakeFiles/planes.dir/flags.make
lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o: ../lib/planes/src/genCutPlanes.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /playpen/software/sreps/Pablo2/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/planes.dir/src/genCutPlanes.cpp.o -c /playpen/software/sreps/Pablo2/lib/planes/src/genCutPlanes.cpp

lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/planes.dir/src/genCutPlanes.cpp.i"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /playpen/software/sreps/Pablo2/lib/planes/src/genCutPlanes.cpp > CMakeFiles/planes.dir/src/genCutPlanes.cpp.i

lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/planes.dir/src/genCutPlanes.cpp.s"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /playpen/software/sreps/Pablo2/lib/planes/src/genCutPlanes.cpp -o CMakeFiles/planes.dir/src/genCutPlanes.cpp.s

lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.requires:
.PHONY : lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.requires

lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.provides: lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.requires
	$(MAKE) -f lib/planes/CMakeFiles/planes.dir/build.make lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.provides.build
.PHONY : lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.provides

lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.provides.build: lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o

# Object files for target planes
planes_OBJECTS = \
"CMakeFiles/planes.dir/src/CutPlanes.cpp.o" \
"CMakeFiles/planes.dir/src/drawAtomCut.cpp.o" \
"CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o" \
"CMakeFiles/planes.dir/src/genCutPlanes.cpp.o"

# External object files for target planes
planes_EXTERNAL_OBJECTS =

libplanes.a: lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o
libplanes.a: lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o
libplanes.a: lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o
libplanes.a: lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o
libplanes.a: lib/planes/CMakeFiles/planes.dir/build.make
libplanes.a: lib/planes/CMakeFiles/planes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../libplanes.a"
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && $(CMAKE_COMMAND) -P CMakeFiles/planes.dir/cmake_clean_target.cmake
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/planes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/planes/CMakeFiles/planes.dir/build: libplanes.a
.PHONY : lib/planes/CMakeFiles/planes.dir/build

lib/planes/CMakeFiles/planes.dir/requires: lib/planes/CMakeFiles/planes.dir/src/CutPlanes.cpp.o.requires
lib/planes/CMakeFiles/planes.dir/requires: lib/planes/CMakeFiles/planes.dir/src/drawAtomCut.cpp.o.requires
lib/planes/CMakeFiles/planes.dir/requires: lib/planes/CMakeFiles/planes.dir/src/P3DCutPlaneView.cpp.o.requires
lib/planes/CMakeFiles/planes.dir/requires: lib/planes/CMakeFiles/planes.dir/src/genCutPlanes.cpp.o.requires
.PHONY : lib/planes/CMakeFiles/planes.dir/requires

lib/planes/CMakeFiles/planes.dir/clean:
	cd /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes && $(CMAKE_COMMAND) -P CMakeFiles/planes.dir/cmake_clean.cmake
.PHONY : lib/planes/CMakeFiles/planes.dir/clean

lib/planes/CMakeFiles/planes.dir/depend:
	cd /playpen/software/sreps/Pablo2/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /playpen/software/sreps/Pablo2 /playpen/software/sreps/Pablo2/lib/planes /playpen/software/sreps/Pablo2/cmake-build-debug /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes /playpen/software/sreps/Pablo2/cmake-build-debug/lib/planes/CMakeFiles/planes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/planes/CMakeFiles/planes.dir/depend

