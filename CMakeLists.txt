project( SOLVERS C CXX )
cmake_minimum_required( VERSION 2.8.8 )

add_library( solverslib 
	src/defs.hpp
	src/emac24.hpp
	src/emac24.cpp
	src/field.hpp
	src/mac22.hpp
	src/mac22.cpp
	src/rk48.hpp
	src/rk48.cpp
	src/solver.hpp )
	
add_executable( solvers src/main.cpp )
target_link_libraries( solvers solverslib )


# TESTING

enable_testing()
find_package( GTest REQUIRED )
include_directories( ${GTEST_INCLUDE_DIRS} )

add_executable( deriv-mac22 src/test/deriv-mac22.cpp )
target_link_libraries( deriv-mac22 solverslib ${GTEST_LIBRARIES} )

add_test( DERIV-F1 deriv-mac22 )
