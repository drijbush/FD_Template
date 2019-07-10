# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(FiniteD Fortran)

add_library(fd constants.f90 fd_template.f90 grid_vectors.f90 laplacian_fd.f90)

if (OpenMP_Fortran_FOUND)
  add_definitions(${OpenMP_Fortran_FLAGS})
endif()

add_executable(test_fd test.f90)
target_link_libraries(test_fd fd m)

find_package(OpenMP)

