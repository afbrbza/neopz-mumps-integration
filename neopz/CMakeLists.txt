cmake_minimum_required(VERSION 3.14)
project(NeoPZMumpsIntegration)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../../neopz_install/ ${CMAKE_SOURCE_DIR}/../../neopz_install/)

add_library(pzmumps STATIC)

target_link_libraries(pzmumps PRIVATE NeoPZ::pz)
target_compile_options(pzmumps PRIVATE -O3 -march=native -fopenmp)
target_link_options(pzmumps PRIVATE -fopenmp)

include(cmake/EnableMUMPS.cmake)
enable_mumps(pzmumps)

add_subdirectory(Solvers)
add_subdirectory(StrMatrix)
add_subdirectory(Matrix)