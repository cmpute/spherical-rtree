cmake_minimum_required(VERSION 3.12)

project(spherical_rtree)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)
include_directories(${NumPy_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_LIST_DIR}/sphertree)

add_subdirectory(sphertree)
