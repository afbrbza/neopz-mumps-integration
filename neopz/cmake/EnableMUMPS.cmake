# ============================================================================
#  EnableMUMPS.cmake
#  Integração do MUMPS SEQUENCIAL (dmumps) ao NeoPZ
#  Baseado no pacote:
#    https://github.com/giavancini/mumps
# ============================================================================

function(enable_mumps target)

    message(STATUS "[EnableMUMPS] Looking for MUMPS (sequential, double)...")

    # evita teste de OpenMP Fortran pelo find_package
    if(NOT DEFINED OpenMP_Fortran_FOUND)
        set(OpenMP_Fortran_FOUND TRUE CACHE INTERNAL "")
        set(OpenMP_Fortran_FLAGS "" CACHE INTERNAL "")
        set(OpenMP_Fortran_LIBRARIES "" CACHE INTERNAL "")
    endif()

    list(APPEND CMAKE_PREFIX_PATH "/opt/mumps")

    # ------------------------------------------------------------------------
    # 1. Localiza o pacote MUMPS
    # ------------------------------------------------------------------------
    find_package(MUMPS CONFIG REQUIRED)

    # ------------------------------------------------------------------------
    # 2. Garante que NÃO estamos usando MUMPS paralelo
    # ------------------------------------------------------------------------
    if(MUMPS_parallel)
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS parallel detected, but this integration "
            "expects SEQUENTIAL MUMPS only."
        )
    endif()

    # ------------------------------------------------------------------------
    # 3. Garante suporte a double precision (dmumps)
    # ------------------------------------------------------------------------
    if(NOT MUMPS_DOUBLE)
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS was not built with DOUBLE precision (dmumps)."
        )
    endif()

    # Força OpenMP para paralelização mesmo com MUMPS sequencial (sem MPI)
    find_package(OpenMP REQUIRED COMPONENTS CXX C)

    message(STATUS "MUMPS_DIR: ${MUMPS_DIR}
        parallel: ${MUMPS_parallel}
        Scotch: ${MUMPS_scotch}
        METIS: ${MUMPS_metis}
        ParMETIS: ${MUMPS_parmetis}
        OpenMP: ${MUMPS_openmp}
    ")

    message(STATUS "MUMPS_intsize64: ${MUMPS_intsize64}
        LAPACK vendor: ${MUMPS_LAPACK_VENDOR}
        Scalapack vendor: ${MUMPS_SCALAPACK_VENDOR}
        real32: ${MUMPS_SINGLE}
        real64: ${MUMPS_DOUBLE}
        complex64: ${MUMPS_COMPLEX}
        complex128: ${MUMPS_COMPLEX16}
    ")

    message(STATUS "[EnableMUMPS] Using SEQUENTIAL dmumps")
    message(STATUS "[EnableMUMPS] MUMPS version: ${MUMPS_VERSION}")

    # ------------------------------------------------------------------------
    # 4. Macros de compilação
    # ------------------------------------------------------------------------
    target_compile_definitions(${target}
        PRIVATE USING_MUMPS
        INTERFACE PZ_USING_MUMPS
    )

    # ------------------------------------------------------------------------
    # 5. Inclui diretórios de headers
    # (o target MUMPS::MUMPS já carrega os includes, mas mantemos explícito)
    # ------------------------------------------------------------------------
    if(DEFINED MUMPS_INCLUDE_DIR)
        target_include_directories(${target}
            PRIVATE ${MUMPS_INCLUDE_DIR}
        )
    endif()

    # ------------------------------------------------------------------------
    # 6. Linka MUMPS com OpenMP (sem MPI paralelo, usa libmpiseq)
    # ------------------------------------------------------------------------
    target_link_libraries(${target}
        PUBLIC MUMPS::MUMPS
        OpenMP::OpenMP_CXX
        OpenMP::OpenMP_C
        "$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
        MPI::MPI_Fortran MPI::MPI_C
    )
    # always link MPI (MPI::MPI_Fortran and MPI::MPI_C) as we use the libmpiseq if not MUMPS_parallel
    
    # Força flags de OpenMP na compilação
    target_compile_options(${target} PRIVATE -fopenmp)
    target_link_options(${target} PRIVATE -fopenmp)

    message(STATUS "[EnableMUMPS] MUMPS sequential enabled for target: ${target}")

endfunction()
