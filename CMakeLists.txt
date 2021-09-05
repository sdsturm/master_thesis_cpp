cmake_minimum_required(VERSION 3.21.0)

project(master_thesis)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Dependencies needed for main library and applications.
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS
    unit_test_framework
    timer
    iostreams	# For gnuplot-iostream.
    system	# For gnuplot-iostream.
    filesystem  # For gnuplot-iostream.
    )

# Add headers of the main library to the include path.
include_directories(./lib/include)

add_subdirectory(./submodules)
add_subdirectory(./lib)

set(GNUPLOT_IOSTREAM_LINK_LIBS Boost::iostreams Boost::system Boost::filesystem)

add_subdirectory(./examples)
add_subdirectory(./production)
