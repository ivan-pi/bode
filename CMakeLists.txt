cmake_minimum_required(VERSION 3.12.0 FATAL_ERROR)

project(bode VERSION 0.1.0 LANGUAGES Fortran C)

add_library(bode 
	src/bode_module.f90
	src/bode_nonlin.f90
	src/bode_solve.f
	src/ltri.f90)

