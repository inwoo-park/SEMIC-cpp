cmake_minimum_required(VERSION 3.0...3.20)
project(SEMIC_cpp VERSION 0.1.0 LANGUAGES CXX Fortran)

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

# enable cmake testing
enable_testing()

# Set C++ standard.
set(CMAKE_CXX_STANDARD 11)

# Initialize variable.
set(LibPySemic libpysemic)
set(LibExample libexample)

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
	# main module for pysemic
	pybind11_add_module(${LibPySemic} src/libpysemic.cpp src/SurfaceEnergyBalance.cpp)
	target_compile_definitions(${LibPySemic} PRIVATE HAVE_PYBIND11)
	if (USE_DEBUG)
		target_compile_definitions(${LibPySemic} PRIVATE HAVE_DEBUG)
	endif()
	if (USE_OPENMP AND OpenMP_CXX_FOUND)
		target_link_libraries(${LibPySemic} PRIVATE OpenMP::OpenMP_CXX)
		target_link_libraries(SemicLib PRIVATE OpenMP::OpenMP_CXX)
	endif()

	# example module for test
	pybind11_add_module(${LibExample} src/libexample.cpp src/libexample.cpp)
	if (USE_OPENMP AND OpenMP_CXX_FOUND)
		target_link_libraries(${LibExample} PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()

# Compile executable SEMIC with
add_executable(SEMIC src/main.cpp)
target_include_directories(SEMIC PUBLIC src)
target_link_libraries(SEMIC SemicLib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include_directories(src)

# Compile original SEMIC written in fortran
add_library(obj_semic_fortran OBJECT
	src/semic-fortran/surface_physics.f90
	src/semic-fortran/utils.f90
)
add_executable(test_semic_fortran_example test/fortran/example.f90
	$<TARGET_OBJECTS:obj_semic_fortran>
)

# Test semic-fortran test.
add_test(NAME TestSemicFortran COMMAND test_semic_fortran_example example.namelist)
set_tests_properties(TestSemicFortran
	PROPERTIES
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/fortran/
)
# Test semic-fortran test: Albedo with Denby
add_test(NAME TestSemicFortranDenby COMMAND test_semic_fortran_example example-denby.namelist)
set_tests_properties(TestSemicFortranDenby
	PROPERTIES
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/fortran/
)

add_test(NAME TestSemicFortranSlater COMMAND test_semic_fortran_example example-slater.namelist)
set_tests_properties(TestSemicFortranSlater
	PROPERTIES
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/fortran/
)

# Test cpp test.
add_executable(test_memory_leakage test/cpp/test_memory_leakage.cpp)
add_test(NAME TestMemoryLeakage COMMAND test_memory_leakage)
target_link_libraries(test_memory_leakage SemicLib)
target_include_directories(test_memory_leakage PUBLIC src)

# Install
# install binary file
#install(TARGETS SEMIC  DESTINATION bin)
# install library file

if (use_pybind11)
	install(TARGETS ${LibPySemic}
		LIBRARY DESTINATION pyseb
	)
	install(TARGETS ${LibExample}
		LIBRARY DESTINATION pyseb
	)
	#install(TARGETS ${LibPySemic} DESTINATION ${CMAKE_SOURCE_DIR}/pyseb/)
	#install(TARGETS ${libexample} DESTINATION ${CMAKE_SOURCE_DIR}/pyseb)
endif()
