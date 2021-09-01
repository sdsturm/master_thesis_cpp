cmake_minimum_required(VERSION 3.21.2)

project(master_thesis)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(Armadillo REQUIRED)

include_directories(./lib/include)

add_subdirectory(./submodules)
add_subdirectory(./lib)
add_subdirectory(./examples)
add_subdirectory(./production)
