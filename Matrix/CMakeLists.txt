# @file neopz/Matrix/CMakeLists.txt  -- CMake file for the Matrix module

target_include_directories(firsttest PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Matrix>
)

set(headers
  TPZSYSMPMumps.h
  TPZYSMPMumps.h
  # TPZSYSMPSpooles.h
  # TPZYSMPSpooles.h
)

set(sources
  TPZSYSMPMumps.cpp
  TPZYSMPMumps.cpp
  # TPZSYSMPSpooles.cpp
  # TPZYSMPSpooles.cpp
)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Matrix)

target_sources(firsttest PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
  list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
