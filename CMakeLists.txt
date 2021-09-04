cmake_minimum_required(VERSION 3.21.0)

project(master_thesis)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(Armadillo REQUIRED)

find_package(Boost REQUIRED COMPONENTS timer iostreams system filesystem)
set(GNUPLOT_IOSTREAM_LINK_LIBS Boost::iostreams Boost::system Boost::filesystem)

include_directories(./lib/include)

add_subdirectory(./submodules)
add_subdirectory(./lib)
add_subdirectory(./examples)
add_subdirectory(./production)
