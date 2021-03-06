
project(VideoNLBayes C CXX Fortran)

message(STATUS "Configuring main project: " ${CMAKE_PROJECT_NAME})

set(PROJECT_VERSION_MAJOR  3)
set(PROJECT_VERSION_MINOR  2)
set(PROJECT_VERSION_PATCH  0)

set(PROJECT_VERSION
    "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
    CACHE INTERNAL "The project version.")

###############################################################################

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")

set(CMAKE_BUILD_TYPE Release)

###############################################################################

cmake_minimum_required(VERSION 2.8.2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
set(CMAKE_VERBOSE_MAKEFILE OFF)

###############################################################################

find_package( CBLAS REQUIRED )

find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTWF_LIBRARIES NAMES fftw3f)

###############################################################################

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

add_subdirectory(lib/iio)
add_subdirectory(src/core)
add_subdirectory(src)

