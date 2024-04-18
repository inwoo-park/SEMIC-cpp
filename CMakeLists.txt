cmake_minimum_required(VERSION 3.0.0)
project(SEMIC-cpp VERSION 0.1.0 LANGUAGES C CXX)

# include(CTest)
# enable_testing()
set(CMAKE_CXX_STANDARD 11)
set(PYSEMIC pySEMIC)

add_subdirectory(pybind11)
# find_package(pybind11 REQUIRED)
pybind11_add_module(${PYSEMIC} src/pySEMIC.cpp src/SurfaceEnergyBalance.cpp)

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
