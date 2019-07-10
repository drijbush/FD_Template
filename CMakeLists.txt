# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(FiniteD Fortran)

add_library(fd_lib constants.f90 FD_template.f90 grid_vectors.f90 laplacian_fd.f90)

add_executable(test_fd test.f90)
target_link_libraries(test_fd fd_lib m)
