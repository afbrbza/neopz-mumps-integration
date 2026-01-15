cmake_minimum_required(VERSION 3.14)
project(NeoPZMumpsIntegration)

set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading layer" FORCE)
find_package(MKL REQUIRED)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../../neopz_install/ ${CMAKE_SOURCE_DIR}/../../neopz_install/)

add_library(pzmumps STATIC)

target_link_libraries(pzmumps PRIVATE NeoPZ::pz MKL::mkl_intel_ilp64 MKL::mkl_intel_thread MKL::mkl_core iomp5)
target_compile_definitions(pzmumps PRIVATE MKL_ILP64)
target_compile_options(pzmumps PRIVATE -O3 -march=native -fopenmp)
target_link_options(pzmumps PRIVATE -fopenmp)

include(cmake/EnableMUMPS.cmake)
enable_mumps(pzmumps)

add_subdirectory(Solvers)
add_subdirectory(StrMatrix)
add_subdirectory(Matrix)