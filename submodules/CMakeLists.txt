# Dependencies of complex_bessel
enable_language(Fortran)    # For submodule complex_bessel
find_package(HDF5 REQUIRED)
find_package(GTest REQUIRED)

add_subdirectory(./complex_bessel)
