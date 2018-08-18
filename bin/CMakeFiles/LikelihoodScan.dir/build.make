# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nicolin/CERN/LimitSetting/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nicolin/CERN/LimitSetting/bin

# Include any dependencies generated for this target.
include CMakeFiles/LikelihoodScan.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LikelihoodScan.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LikelihoodScan.dir/flags.make

CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.o: CMakeFiles/LikelihoodScan.dir/flags.make
CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.o: /home/nicolin/CERN/LimitSetting/src/LikelihoodScan.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolin/CERN/LimitSetting/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.o -c /home/nicolin/CERN/LimitSetting/src/LikelihoodScan.C

CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolin/CERN/LimitSetting/src/LikelihoodScan.C > CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.i

CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolin/CERN/LimitSetting/src/LikelihoodScan.C -o CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.s

# Object files for target LikelihoodScan
LikelihoodScan_OBJECTS = \
"CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.o"

# External object files for target LikelihoodScan
LikelihoodScan_EXTERNAL_OBJECTS =

LikelihoodScan: CMakeFiles/LikelihoodScan.dir/LikelihoodScan.C.o
LikelihoodScan: CMakeFiles/LikelihoodScan.dir/build.make
LikelihoodScan: libWorkspaceCalculator.a
LikelihoodScan: libProfileLikelihoodTestStatEnhanced.a
LikelihoodScan: libLikelihoodCalculator.a
LikelihoodScan: /usr/local/lib/libCore.so
LikelihoodScan: /usr/local/lib/libImt.so
LikelihoodScan: /usr/local/lib/libRIO.so
LikelihoodScan: /usr/local/lib/libNet.so
LikelihoodScan: /usr/local/lib/libHist.so
LikelihoodScan: /usr/local/lib/libGraf.so
LikelihoodScan: /usr/local/lib/libGraf3d.so
LikelihoodScan: /usr/local/lib/libGpad.so
LikelihoodScan: /usr/local/lib/libTree.so
LikelihoodScan: /usr/local/lib/libTreePlayer.so
LikelihoodScan: /usr/local/lib/libRint.so
LikelihoodScan: /usr/local/lib/libPostscript.so
LikelihoodScan: /usr/local/lib/libMatrix.so
LikelihoodScan: /usr/local/lib/libPhysics.so
LikelihoodScan: /usr/local/lib/libMathCore.so
LikelihoodScan: /usr/local/lib/libThread.so
LikelihoodScan: /usr/local/lib/libMultiProc.so
LikelihoodScan: /usr/local/lib/libRooStats.so
LikelihoodScan: /usr/local/lib/libRooFit.so
LikelihoodScan: /usr/local/lib/libRooFitCore.so
LikelihoodScan: /usr/local/lib/libHistFactory.so
LikelihoodScan: /usr/local/lib/libMinuit2.so
LikelihoodScan: /usr/local/lib/libHist.so
LikelihoodScan: /usr/local/lib/libMatrix.so
LikelihoodScan: /usr/local/lib/libMathCore.so
LikelihoodScan: /usr/local/lib/libImt.so
LikelihoodScan: /usr/local/lib/libRIO.so
LikelihoodScan: /usr/local/lib/libThread.so
LikelihoodScan: /usr/local/lib/libCore.so
LikelihoodScan: CMakeFiles/LikelihoodScan.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nicolin/CERN/LimitSetting/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable LikelihoodScan"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LikelihoodScan.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LikelihoodScan.dir/build: LikelihoodScan

.PHONY : CMakeFiles/LikelihoodScan.dir/build

CMakeFiles/LikelihoodScan.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LikelihoodScan.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LikelihoodScan.dir/clean

CMakeFiles/LikelihoodScan.dir/depend:
	cd /home/nicolin/CERN/LimitSetting/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nicolin/CERN/LimitSetting/src /home/nicolin/CERN/LimitSetting/src /home/nicolin/CERN/LimitSetting/bin /home/nicolin/CERN/LimitSetting/bin /home/nicolin/CERN/LimitSetting/bin/CMakeFiles/LikelihoodScan.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LikelihoodScan.dir/depend
