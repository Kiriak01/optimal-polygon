# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/george/Desktop/optimal-algo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/george/Desktop/optimal-algo

# Include any dependencies generated for this target.
include includes/CMakeFiles/opt_algorithms.dir/depend.make

# Include the progress variables for this target.
include includes/CMakeFiles/opt_algorithms.dir/progress.make

# Include the compile flags for this target's objects.
include includes/CMakeFiles/opt_algorithms.dir/flags.make

includes/CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.o: includes/CMakeFiles/opt_algorithms.dir/flags.make
includes/CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.o: includes/opt_algorithms.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/george/Desktop/optimal-algo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object includes/CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.o"
	cd /home/george/Desktop/optimal-algo/includes && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.o -c /home/george/Desktop/optimal-algo/includes/opt_algorithms.cpp

includes/CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.i"
	cd /home/george/Desktop/optimal-algo/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/george/Desktop/optimal-algo/includes/opt_algorithms.cpp > CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.i

includes/CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.s"
	cd /home/george/Desktop/optimal-algo/includes && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/george/Desktop/optimal-algo/includes/opt_algorithms.cpp -o CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.s

# Object files for target opt_algorithms
opt_algorithms_OBJECTS = \
"CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.o"

# External object files for target opt_algorithms
opt_algorithms_EXTERNAL_OBJECTS =

includes/libopt_algorithms.a: includes/CMakeFiles/opt_algorithms.dir/opt_algorithms.cpp.o
includes/libopt_algorithms.a: includes/CMakeFiles/opt_algorithms.dir/build.make
includes/libopt_algorithms.a: includes/CMakeFiles/opt_algorithms.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/george/Desktop/optimal-algo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libopt_algorithms.a"
	cd /home/george/Desktop/optimal-algo/includes && $(CMAKE_COMMAND) -P CMakeFiles/opt_algorithms.dir/cmake_clean_target.cmake
	cd /home/george/Desktop/optimal-algo/includes && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/opt_algorithms.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
includes/CMakeFiles/opt_algorithms.dir/build: includes/libopt_algorithms.a

.PHONY : includes/CMakeFiles/opt_algorithms.dir/build

includes/CMakeFiles/opt_algorithms.dir/clean:
	cd /home/george/Desktop/optimal-algo/includes && $(CMAKE_COMMAND) -P CMakeFiles/opt_algorithms.dir/cmake_clean.cmake
.PHONY : includes/CMakeFiles/opt_algorithms.dir/clean

includes/CMakeFiles/opt_algorithms.dir/depend:
	cd /home/george/Desktop/optimal-algo && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/george/Desktop/optimal-algo /home/george/Desktop/optimal-algo/includes /home/george/Desktop/optimal-algo /home/george/Desktop/optimal-algo/includes /home/george/Desktop/optimal-algo/includes/CMakeFiles/opt_algorithms.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : includes/CMakeFiles/opt_algorithms.dir/depend

