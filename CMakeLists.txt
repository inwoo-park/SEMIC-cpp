cmake_minimum_required(VERSION 3.0...3.20)
project(SEMIC_cpp VERSION 0.1.0 LANGUAGES C CXX)

# Include default cmake package.
include(CMakePrintHelpers)

# Set compiler dependig on system
# if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# 	#set(CMAKE_CXX_COMPILER "C:/")
# elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
# 	set(CMAKE_C_COMPILE "/usr/bin/gcc")
# 	set(CMAKE_CXX_COMPILER "/usr/bin/g++")
# else()
# 	message(FATAL_ERROR "Unsupported operating system.")
# endif()

# Set C++ standard.
set(CMAKE_CXX_STANDARD 11)

# Initialize variable.
set(PYSEMIC libpysemic)

option(USE_OPENMP "Use OpenMP for parallelization" ON)
option(use_pybind11 "Use pybind11" ON)
option(USE_DEBUG "Use debugging" OFF)

# Find OpenMP package.
if (USE_OPENMP)
	find_package(OpenMP REQUIRED)
	if (NOT OpenMP_CXX_FOUND)
		message(WARNING "OpenMP not found, proceeding without it.")
	endif()
endif()

# Find pybind11 package.
if (use_pybind11)
	if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
		find_package(Python3 REQUIRED Interpreter Development)
		cmake_print_variables(Python3_EXECUTABLE)

		execute_process(
		COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
			OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
		cmake_print_variables(PYBIND11_CMAKE_DIR)
		find_package(pybind11 REQUIRED PATHS ${PYBIND11_CMAKE_DIR})
	elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
		find_package(Python3 REQUIRED Interpreter Development)
		cmake_print_variables(Python3_EXECUTABLE)

		execute_process(
		COMMAND ${Python3_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
			OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
		cmake_print_variables(PYBIND11_CMAKE_DIR)
		find_package(pybind11 REQUIRED PATHS ${PYBIND11_CMAKE_DIR})
	endif()
endif()

# Initialize SEMIC library for main file
add_library(SemicLib STATIC src/SurfaceEnergyBalance.cpp
 src/SurfaceEnergyBalance.h
 src/SemicParameters.h
 src/DataArray.h
 )

if (use_pybind11)
	pybind11_add_module(${PYSEMIC} src/libpysemic.cpp src/SurfaceEnergyBalance.cpp)
	target_compile_definitions(${PYSEMIC} PRIVATE HAVE_PYBIND11)
	if (USE_DEBUG)
		target_compile_definitions(${PYSEMIC} PRIVATE HAVE_DEBUG)
	endif()
	if (USE_OPENMP AND OpenMP_CXX_FOUND)
		target_link_libraries(${PYSEMIC} PRIVATE OpenMP::OpenMP_CXX)
		target_link_libraries(SemicLib PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()

# Initialize SEMIC
add_executable(SEMIC src/main.cpp)
target_include_directories(SEMIC PUBLIC src)
target_link_libraries(SEMIC SemicLib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include_directories(src)

# Test cpp test.
add_executable(test_memory_leakage test/cpp/test_memory_leakage.cpp)
target_link_libraries(test_memory_leakage SemicLib)
target_include_directories(test_memory_leakage PUBLIC src)

# Install
# install binary file
install(TARGETS SEMIC  DESTINATION bin)
# install library file

if (use_pybind11)
	install(TARGETS ${PYSEMIC} DESTINATION lib)
	install(TARGETS ${PYSEMIC} DESTINATION ../pyseb/)
endif()
