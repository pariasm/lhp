cmake_policy(SET CMP0048 NEW)

project(vnlmeans C)

message(STATUS "Configuring main project: " ${CMAKE_PROJECT_NAME})

###############################################################################

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -ggdb -fsanitize=address")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -ggdb")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall -ggdb -O6 -w -s -ffast-math -fno-finite-math-only -fomit-frame-pointer -fstrength-reduce")

set(CMAKE_BUILD_TYPE Release)

###############################################################################

cmake_minimum_required(VERSION 2.8.2)
set(CMAKE_VERBOSE_MAKEFILE OFF)

###############################################################################

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

add_subdirectory(lib/argparse) # command line parsing
add_subdirectory(lib/iio)      # image i/o
add_subdirectory(lib/tvl1flow) # tvl1 optical flow
add_subdirectory(src)
