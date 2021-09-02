cmake_minimum_required(VERSION 3.21.2)

project(master_thesis)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
find_package(Armadillo REQUIRED)

set(GNUPLOT_IOSTREAM_LINK_LIBS Boost::iostreams Boost::system Boost::filesystem)

include_directories(./lib/include)

# Directory for data files to use in TikZ/PGFPlots.
add_compile_definitions(MTHESIS_RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/results/")

add_subdirectory(./submodules)
add_subdirectory(./lib)
add_subdirectory(./examples)
add_subdirectory(./production)
