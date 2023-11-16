# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_24_1/Linux64bit+3.10-2.17/bin/cmake

# The command to remove a file.
RM = /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_24_1/Linux64bit+3.10-2.17/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff

# Include any dependencies generated for this target.
include src/CMakeFiles/geoEff.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/geoEff.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/geoEff.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/geoEff.dir/flags.make

src/CMakeFiles/geoEff.dir/geoEff.cpp.o: src/CMakeFiles/geoEff.dir/flags.make
src/CMakeFiles/geoEff.dir/geoEff.cpp.o: src/geoEff.cpp
src/CMakeFiles/geoEff.dir/geoEff.cpp.o: src/CMakeFiles/geoEff.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/geoEff.dir/geoEff.cpp.o"
	cd /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src && /cvmfs/larsoft.opensciencegrid.org/products/gcc/v9_3_0/Linux64bit+3.10-2.17/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/geoEff.dir/geoEff.cpp.o -MF CMakeFiles/geoEff.dir/geoEff.cpp.o.d -o CMakeFiles/geoEff.dir/geoEff.cpp.o -c /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src/geoEff.cpp

src/CMakeFiles/geoEff.dir/geoEff.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/geoEff.dir/geoEff.cpp.i"
	cd /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src && /cvmfs/larsoft.opensciencegrid.org/products/gcc/v9_3_0/Linux64bit+3.10-2.17/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src/geoEff.cpp > CMakeFiles/geoEff.dir/geoEff.cpp.i

src/CMakeFiles/geoEff.dir/geoEff.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/geoEff.dir/geoEff.cpp.s"
	cd /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src && /cvmfs/larsoft.opensciencegrid.org/products/gcc/v9_3_0/Linux64bit+3.10-2.17/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src/geoEff.cpp -o CMakeFiles/geoEff.dir/geoEff.cpp.s

# Object files for target geoEff
geoEff_OBJECTS = \
"CMakeFiles/geoEff.dir/geoEff.cpp.o"

# External object files for target geoEff
geoEff_EXTERNAL_OBJECTS =

lib/libgeoEff.so: src/CMakeFiles/geoEff.dir/geoEff.cpp.o
lib/libgeoEff.so: src/CMakeFiles/geoEff.dir/build.make
lib/libgeoEff.so: src/CMakeFiles/geoEff.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../lib/libgeoEff.so"
	cd /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/geoEff.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/geoEff.dir/build: lib/libgeoEff.so
.PHONY : src/CMakeFiles/geoEff.dir/build

src/CMakeFiles/geoEff.dir/clean:
	cd /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src && $(CMAKE_COMMAND) -P CMakeFiles/geoEff.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/geoEff.dir/clean

src/CMakeFiles/geoEff.dir/depend:
	cd /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src /dune/app/users/awilkins/nd-sim-tools/inputs/DUNE_ND_GeoEff/src/CMakeFiles/geoEff.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/geoEff.dir/depend

