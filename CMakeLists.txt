cmake_minimum_required(VERSION 3.21.0)

project(master-thesis)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Dependencies at top level.
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS
    timer	# For performance analysis.
    iostreams	# For gnuplot-iostream.
    system	# For gnuplot-iostream.
    filesystem  # For gnuplot-iostream.
    )

# Libraries needed by gnuplot-iostream
set(GNUPLOT_IOSTREAM_LINK_LIBS Boost::iostreams Boost::system Boost::filesystem)


# Structure of this project.
add_subdirectory(./submodules)	 	# External modules.

include_directories(./lib/include) 	# Main library include directory.
add_subdirectory(./lib)			# Main library folder.

add_subdirectory(./examples) 		# Illustrative examples.
add_subdirectory(./production)		# The real applications.
