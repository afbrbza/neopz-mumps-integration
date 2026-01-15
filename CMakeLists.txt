cmake_minimum_required(VERSION 3.14)
project(FirstTest)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/../neopz_install/)

set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading layer" FORCE)
find_package(MKL REQUIRED)

################################################################
# PARDISO example
add_executable(main.pardiso main.pardiso.cpp)
# MUMPS example
add_executable(main.mumps main.mumps.cpp)
################################################################
target_sources(main.pardiso PRIVATE ${CMAKE_SOURCE_DIR}/OutputCapture.cpp)
target_sources(main.mumps PRIVATE ${CMAKE_SOURCE_DIR}/OutputCapture.cpp)

add_subdirectory(neopz)

target_link_libraries(main.pardiso PRIVATE NeoPZ::pz MKL::mkl_intel_ilp64 MKL::mkl_intel_thread MKL::mkl_core iomp5)
target_compile_definitions(main.pardiso PRIVATE MKL_ILP64)
target_compile_options(main.pardiso PRIVATE -O3 -march=native -fopenmp)
target_link_options(main.pardiso PRIVATE -fopenmp)

target_link_libraries(main.mumps PRIVATE NeoPZ::pz pzmumps)