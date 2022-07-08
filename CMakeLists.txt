cmake_minimum_required(VERSION 3.21.0)

project(master-thesis)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS timer)

add_subdirectory(./submodules)	 	# External modules.

include_directories(./lib/include) 	# Main library include directory.
add_subdirectory(./lib)			# Main library folder.

add_subdirectory(./examples) 		# Illustrative examples.
add_subdirectory(./production)		# The real applications.
