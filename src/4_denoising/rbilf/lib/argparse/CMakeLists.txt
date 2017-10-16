cmake_minimum_required(VERSION 2.8)

project(argparse C)

message(STATUS "Configuring project: " ${CMAKE_PROJECT_NAME})

###############################################################################

# NOTE: -fPIC is added so that the static library can be included into
#       shared libraries.
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -fPIC -O3 -g -ggdb -Wall -Wextra")
set(CMAKE_BUILD_TYPE Release)

###############################################################################

add_library(argparse STATIC argparse.c)
add_library(argparse_shared SHARED argparse.c)

target_link_libraries(argparse m)
target_link_libraries(argparse_shared m)

###############################################################################

# Expose public includes to other subprojects through cache variable
set(ARGPARSE_INCLUDE_DIR ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

