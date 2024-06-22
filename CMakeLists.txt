cmake_minimum_required(VERSION 3.0...3.28)
project(SEMIC-cpp VERSION 0.1.0 LANGUAGES C CXX)

# Set compiler dependig on system
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	#set(CMAKE_CXX_COMPILER "C:/")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	set(CMAKE_C_COMPILE "/usr/bin/gcc")
	set(CMAKE_CXX_COMPILER "/usr/bin/g++")
else()
	message(FATAL_ERROR "Unsupported operating system.")
endif()

# include(CTest)
# enable_testing()
set(CMAKE_CXX_STANDARD 11)
set(PYSEMIC libpysemic)

option(USE_OPENMP "Use OpenMP for parallelization" OFF)

# find openmp package
if (USE_OPENMP)
	find_package(OpenMP)
	if (NOT OpenMP_CXX_FOUND)
		message(WARNING "OpenMP not found, proceeding without it.")
	endif()
endif()

# pybind11 specials
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	add_subdirectory(pybind11)
	#find_package(pybind11 REQUIRED)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	find_package(Python3 REQUIRED Interpreter Development)
	find_package(pybind11 REQUIRED)
endif()

pybind11_add_module(${PYSEMIC} src/pySEMIC.cpp src/SurfaceEnergyBalance.cpp)
if (USE_OPENMP AND OpenMP_CXX_FOUND)
	target_link_libraries(${PYSEMIC} PRIVATE OpenMP::OpenMP_CXX)
endif()

# pybind11_add_module(MyLib src/pybind.cpp)
#target_compile_definitions(pySEMIC PRIVATE)

# Initialize SEMIC library for main file
add_library(SemicLib STATIC src/SurfaceEnergyBalance.cpp
 src/SurfaceEnergyBalance.h
 src/SemicParameters.h)

# Initialize SEMIC
add_executable(SEMIC src/main.cpp)
target_include_directories(SEMIC PUBLIC src)
target_link_libraries(SEMIC SemicLib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include_directories(src)

# Install
# install binary file
install(TARGETS SEMIC  DESTINATION bin)
# install library file
install(TARGETS ${PYSEMIC} DESTINATION lib)
