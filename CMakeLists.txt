cmake_minimum_required(VERSION 3.0.0)
project(SEMIC-cpp VERSION 0.1.0 LANGUAGES C CXX)

# include(CTest)
enable_testing()

add_library(SemicLib STATIC src/SurfaceEnergyBalance.cpp src/SurfaceEnergyBalance.h)
add_executable(SEMIC src/main.cpp)
target_link_libraries(SEMIC SemicLib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include_directories(src)
