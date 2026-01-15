cmake_minimum_required(VERSION 3.14)
project(FirstTest)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/../neopz_install/)

################################################################
# PARDISO example
# add_executable(firsttest main.pardiso.cpp)
# MUMPS example
add_executable(firsttest main.mumps.cpp)
################################################################
target_sources(firsttest PRIVATE ${CMAKE_SOURCE_DIR}/OutputCapture.cpp)

add_subdirectory(neopz)

target_link_libraries(firsttest PRIVATE NeoPZ::pz pzmumps)