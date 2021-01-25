cmake_minimum_required(VERSION 3.3)
project(myspice)

option(EXTENDED OFF)

set(CMAKE_CXX_FLAGS "--std=c++20 -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-builtin -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math")

add_executable(
	myspice
	"${CMAKE_SOURCE_DIR}/src/myspice.cpp"
	"${CMAKE_SOURCE_DIR}/src/mna.cpp"
	"${CMAKE_SOURCE_DIR}/src/circuit.cpp"
	"${CMAKE_SOURCE_DIR}/src/legacy.cpp"
	"${CMAKE_SOURCE_DIR}/src/extended.cpp"
)

if(EXTENDED)
	add_definitions(-DEXTENDED_MODE)
endif()