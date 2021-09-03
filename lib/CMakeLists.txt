project(mthesis-lib)

add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-pedantic)
# add_compile_options(-Werror)

add_library(${PROJECT_NAME}
    ./src/solution_domain.cpp
    ./src/tlgf.cpp
    ./src/spectral_gf.cpp
    ./src/partition_extrapolation.cpp
    ./src/sommerfeld_integrals.cpp
    ./src/sgf.cpp
    ./src/dgf.cpp
    ./src/gpof.cpp
    ./src/dcim_utils.cpp
    ./src/dcim_three_level_v3.cpp
    ./src/fmm.cpp
    )

find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
find_package(Armadillo REQUIRED)

set(BLA_VENDOR OpenBLAS) # Repacement for armadillo runtime library.
find_package(BLAS)

target_compile_definitions(${PROJECT_NAME} PUBLIC ARMA_DONT_USE_WRAPPER)

target_link_libraries(${PROJECT_NAME}
    openblas    # Instead of armadillo runtime library.
    # armadillo
    GSL::gsl GSL::gslcblas
    complex_bessel
    )

add_subdirectory(./test)
