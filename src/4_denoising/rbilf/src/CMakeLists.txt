
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# internal includes
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${ARGPARSE_INCLUDE_DIR})
include_directories(${IIO_INCLUDE_DIR})

add_executable(vnlmeans main.c)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(vnlmeans iio argparse)

