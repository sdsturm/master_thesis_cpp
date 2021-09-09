project(mthesis-lib)

################################################################################
#                                 COMPILATION
################################################################################

add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-pedantic)
# add_compile_options(-Werror)

add_library(${PROJECT_NAME} SHARED
    ./src/solution_domain.cpp
    ./src/tlgf.cpp

    ./src/partition_extrapolation.cpp
    ./src/sommerfeld_integrals.cpp

    ./src/gf/scalar.cpp
    ./src/gf/dyadic.cpp

    ./src/dcim/gpof.cpp
    ./src/dcim/dcim_base.cpp
    ./src/dcim/three_level_v1.cpp
    ./src/dcim/three_level_v2.cpp

    ./src/fmm/helpers.cpp
    ./src/fmm/free_space.cpp
    ./src/fmm/gegenbauer.cpp

    ./src/mom.cpp
    )

################################################################################
#                                   LINKING
################################################################################

# Indispensible libraries.
find_package(Boost REQUIRED) # mthesis-lib usis only header-only parts.
find_package(GSL REQUIRED)
find_package(Armadillo REQUIRED)

target_link_libraries(${PROJECT_NAME}
    GSL::gsl GSL::gslcblas
    complex_bessel
    )

# Link against OpenMP if available.
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

# Linear algebra packages used by Armadillo.
set(BLA_VENDOR OpenBLAS) # Repacement for armadillo runtime library.
find_package(BLAS)
if (BLAS_FOUND) # If openBLAS is available, link armadillo against that.
    target_compile_definitions(${PROJECT_NAME} PUBLIC ARMA_DONT_USE_WRAPPER)
    target_link_libraries(${PROJECT_NAME} openblas)
else() # Otherwise link against armadillo runtime library.
    target_link_libraries(${PROJECT_NAME} armadillo)
endif()

################################################################################
#                                 UNIT TESTS
################################################################################

add_subdirectory(./test)
