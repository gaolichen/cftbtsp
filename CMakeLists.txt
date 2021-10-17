# CMakeLists files in this project can

# refer to the root source directory of the project as ${BTSP_SOURCE_DIR} and

# to the root binary directory of the project as ${BTSP_BINARY_DIR}.

cmake_minimum_required (VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

project (BTSP)
SET (CMAKE_BUILD_TYPE Release)
# SET (CMAKE_BUILD_TYPE Debug)

add_subdirectory (src)

