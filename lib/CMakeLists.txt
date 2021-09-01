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
    ./src/gpof.cpp
    )
target_link_libraries(${PROJECT_NAME}
    armadillo
    complex_bessel
    )

add_subdirectory(./test)
