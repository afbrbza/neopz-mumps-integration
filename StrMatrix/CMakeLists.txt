# @file neopz/StrMatrix/CMakeLists.txt  -- CMake file for the StrMatrix module

target_include_directories(firsttest PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix>
                           )
#@TODO: check if these public headers are exaggerated
set(public_headers
    #sparse
    TPZSSpStructMatrixMumps.h
    TPZSpStructMatrixMumps.h
    TPZSSpStructMatrixSpooles.h
    TPZSpStructMatrixSpooles.h
   )
set(headers
    #sparse
    TPZSSpStructMatrixMumps.h
    TPZSpStructMatrixMumps.h
    TPZSSpStructMatrixSpooles.h
    TPZSpStructMatrixSpooles.h
   )

set(sources
    #sparse
    TPZSSpStructMatrixMumps.cpp
    TPZSpStructMatrixMumps.cpp
    TPZSSpStructMatrixSpooles.cpp
    TPZSpStructMatrixSpooles.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix)

target_sources(firsttest PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)