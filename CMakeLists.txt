cmake_minimum_required(VERSION 3.14)
project(FirstTest)
set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading layer" FORCE)
find_package(MKL REQUIRED)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/../neopz_install/)

add_executable(firsttest main.cpp)
# add_executable(firsttest main.pardiso.cpp)
#############################################################################################################
# target_link_libraries(firsttest PRIVATE NeoPZ::pz)
target_link_libraries(firsttest PRIVATE NeoPZ::pz MKL::mkl_intel_ilp64 MKL::mkl_intel_thread MKL::mkl_core iomp5)
target_compile_definitions(firsttest PRIVATE MKL_ILP64)
target_compile_options(firsttest PRIVATE -O3 -march=native -fopenmp)
target_link_options(firsttest PRIVATE -fopenmp)
#############################################################################################################
# target_link_libraries(firsttest PRIVATE NeoPZ::pz -fsanitize=address)
#############################################################################################################
# add_compile_options(-fsanitize=address -g)

add_subdirectory(Solvers)
add_subdirectory(StrMatrix)
add_subdirectory(Matrix)

#enabling MUMPS library
# option(USING_MUMPS "Whether the MUMPS library will be linked in" OFF)
# if(USING_MUMPS)
include(cmake/EnableMUMPS.cmake)
enable_mumps(firsttest)
# endif()