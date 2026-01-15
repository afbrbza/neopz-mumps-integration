# ============================================================================
#  EnableSPOOLES.cmake
#  Integração do SPOOLES ao projeto FirstTest
# ============================================================================

function(enable_spooles target)

    message(STATUS "[EnableSPOOLES] Looking for SPOOLES library...")

    # ------------------------------------------------------------------------
    # 1. Define SPOOLES installation path (pode ser ajustado pelo usuário)
    # ------------------------------------------------------------------------
    if(NOT DEFINED SPOOLES_DIR)
        set(SPOOLES_DIR "/opt/spooles" CACHE PATH "SPOOLES installation directory")
    endif()

    # Alternate search paths
    set(SPOOLES_SEARCH_PATHS
        ${SPOOLES_DIR}
        /opt/spooles
        /usr/local/spooles
        /tmp/spooles
    )

    # ------------------------------------------------------------------------
    # 2. Procura pela biblioteca spooles.a
    # ------------------------------------------------------------------------
    find_library(SPOOLES_LIBRARY
        NAMES spooles.a libspooles.a
        PATHS ${SPOOLES_SEARCH_PATHS}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH
    )

    if(NOT SPOOLES_LIBRARY)
        find_library(SPOOLES_LIBRARY
            NAMES spooles.a libspooles.a
        )
    endif()

    if(NOT SPOOLES_LIBRARY)
        message(FATAL_ERROR
            "[EnableSPOOLES] SPOOLES library not found!\n"
            "Please install SPOOLES using install.spooles.sh script or set SPOOLES_DIR."
        )
    endif()

    message(STATUS "[EnableSPOOLES] Found SPOOLES library: ${SPOOLES_LIBRARY}")

    # ------------------------------------------------------------------------
    # 3. Procura pelos headers do SPOOLES
    # ------------------------------------------------------------------------
    find_path(SPOOLES_INCLUDE_DIR
        NAMES InpMtx.h
        PATHS ${SPOOLES_SEARCH_PATHS}
        PATH_SUFFIXES include . ""
        NO_DEFAULT_PATH
    )

    if(NOT SPOOLES_INCLUDE_DIR)
        # Try to find in the same directory as the library
        get_filename_component(SPOOLES_LIB_DIR ${SPOOLES_LIBRARY} DIRECTORY)
        set(SPOOLES_INCLUDE_DIR ${SPOOLES_LIB_DIR})
    endif()

    if(NOT EXISTS "${SPOOLES_INCLUDE_DIR}/InpMtx.h")
        message(FATAL_ERROR
            "[EnableSPOOLES] SPOOLES headers not found!\n"
            "Expected InpMtx.h at: ${SPOOLES_INCLUDE_DIR}"
        )
    endif()

    message(STATUS "[EnableSPOOLES] Found SPOOLES headers: ${SPOOLES_INCLUDE_DIR}")

    # ------------------------------------------------------------------------
    # 4. Verifica suporte a threads (pthread)
    # ------------------------------------------------------------------------
    find_package(Threads REQUIRED)

    # ------------------------------------------------------------------------
    # 5. Macros de compilação
    # ------------------------------------------------------------------------
    target_compile_definitions(${target}
        PRIVATE USING_SPOOLES
        INTERFACE PZ_USING_SPOOLES
    )

    # Verifica se SPOOLES foi compilado com suporte a threads
    # (procura por MT no diretório)
    if(EXISTS "${SPOOLES_LIB_DIR}/MT")
        message(STATUS "[EnableSPOOLES] SPOOLES thread support detected (MT)")
        target_compile_definitions(${target}
            PRIVATE SPOOLES_THREAD_SUPPORT
        )
    endif()

    # ------------------------------------------------------------------------
    # 6. Inclui diretórios de headers do SPOOLES
    # ------------------------------------------------------------------------
    target_include_directories(${target}
        PRIVATE 
            ${SPOOLES_INCLUDE_DIR}
            ${SPOOLES_INCLUDE_DIR}/A2
            ${SPOOLES_INCLUDE_DIR}/BKL
            ${SPOOLES_INCLUDE_DIR}/BPG
            ${SPOOLES_INCLUDE_DIR}/Chv
            ${SPOOLES_INCLUDE_DIR}/ChvList
            ${SPOOLES_INCLUDE_DIR}/ChvManager
            ${SPOOLES_INCLUDE_DIR}/Coords
            ${SPOOLES_INCLUDE_DIR}/DenseMtx
            ${SPOOLES_INCLUDE_DIR}/DSTree
            ${SPOOLES_INCLUDE_DIR}/Drand
            ${SPOOLES_INCLUDE_DIR}/DV
            ${SPOOLES_INCLUDE_DIR}/EGraph
            ${SPOOLES_INCLUDE_DIR}/ETree
            ${SPOOLES_INCLUDE_DIR}/FrontMtx
            ${SPOOLES_INCLUDE_DIR}/GPart
            ${SPOOLES_INCLUDE_DIR}/Graph
            ${SPOOLES_INCLUDE_DIR}/I2Ohash
            ${SPOOLES_INCLUDE_DIR}/IIheap
            ${SPOOLES_INCLUDE_DIR}/InpMtx
            ${SPOOLES_INCLUDE_DIR}/IV
            ${SPOOLES_INCLUDE_DIR}/IVL
            ${SPOOLES_INCLUDE_DIR}/Ideq
            ${SPOOLES_INCLUDE_DIR}/Lock
            ${SPOOLES_INCLUDE_DIR}/MSMD
            ${SPOOLES_INCLUDE_DIR}/MPI
            ${SPOOLES_INCLUDE_DIR}/MT
            ${SPOOLES_INCLUDE_DIR}/Perm
            ${SPOOLES_INCLUDE_DIR}/QRser
            ${SPOOLES_INCLUDE_DIR}/SemiImplMtx
            ${SPOOLES_INCLUDE_DIR}/SolveMap
            ${SPOOLES_INCLUDE_DIR}/SubMtx
            ${SPOOLES_INCLUDE_DIR}/SubMtxList
            ${SPOOLES_INCLUDE_DIR}/SubMtxManager
            ${SPOOLES_INCLUDE_DIR}/SymbFac
            ${SPOOLES_INCLUDE_DIR}/Tree
            ${SPOOLES_INCLUDE_DIR}/Utilities
            ${SPOOLES_INCLUDE_DIR}/misc
    )

    # ------------------------------------------------------------------------
    # 7. Linka SPOOLES com pthread
    # ------------------------------------------------------------------------
    target_link_libraries(${target}
        PRIVATE 
            ${SPOOLES_LIBRARY}
            Threads::Threads
            m  # math library
    )

    # Força flags de OpenMP na compilação (para paralelização)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
        message(STATUS "[EnableSPOOLES] OpenMP enabled for parallel assembly")
    endif()

    message(STATUS "[EnableSPOOLES] SPOOLES enabled for target: ${target}")

endfunction()
