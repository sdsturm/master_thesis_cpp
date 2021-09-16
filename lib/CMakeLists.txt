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
    ./src/nonspectral.cpp
    ./src/sommerfeld_integrals.cpp

    ./src/gf/scalar.cpp
    ./src/gf/dyadic.cpp

    ./src/dcim/gpof.cpp
    ./src/dcim/dcim_base.cpp
    ./src/dcim/two_level.cpp
    ./src/dcim/three_level_v1.cpp
    ./src/dcim/three_level_v2.cpp
    ./src/dcim/three_level_v3.cpp

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
find_package(Armadillo REQUIRED) # See documentation for dependencies.

target_link_libraries(${PROJECT_NAME}
    armadillo
    GSL::gsl GSL::gslcblas
    complex_bessel
    )

# Link against OpenMP if available.
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

################################################################################
#                                 UNIT TESTS
################################################################################

add_subdirectory(./test)
