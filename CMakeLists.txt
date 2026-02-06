cmake_minimum_required(VERSION 3.14)
project(FirstTest)

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/neopz_install/ ${CMAKE_SOURCE_DIR}/../neopz_install/)

# add_executable(main_mumps main_mumps.cpp OutputCapture.cpp)
# target_link_libraries(main_mumps PRIVATE NeoPZ::pz)

# add_executable(main_tosco main_tosco.cpp)
# target_link_libraries(main_tosco PRIVATE NeoPZ::pz)

add_executable(main_mumps_3d main_mumps_3d.cpp OutputCapture.cpp)
target_link_libraries(main_mumps_3d PRIVATE NeoPZ::pz)
