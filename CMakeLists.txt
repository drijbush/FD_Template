# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(FiniteD Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(FiniteD_FORD_DIR   "${CMAKE_CURRENT_SOURCE_DIR}/ford")
set(FiniteD_DOC_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/doc")
set(FiniteD_SRC_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Add OMP
find_package(OpenMP)
if (OpenMP_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS})
endif()

add_subdirectory(${FiniteD_SRC_DIR})
add_subdirectory(${FiniteD_FORD_DIR})
