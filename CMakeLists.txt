cmake_minimum_required(VERSION 3.12.0 FATAL_ERROR)

project(bode VERSION 0.1.0 LANGUAGES Fortran C)

add_library(bode
	src/bode_module.f90
	src/bode_nonlin.f90
	src/bode_solve.f90)

add_executable(test_banded1 src/test_banded1.f90)
target_link_libraries(test_banded1 bode)

add_executable(test_banded2 src/test_banded2.f90)
target_link_libraries(test_banded2 bode)

add_executable(test src/test.f90)
target_link_libraries(test bode)

add_executable(test_old src/test.f90 src/bode.f)

add_executable(test_vdp src/test_van_der_pol.f90)
target_link_libraries(test_vdp bode)

add_executable(test_vdp_old src/test_van_der_pol.f90 src/bode.f)

# Banded procedures by Jeff Thorson
add_library(band
    src/band.f
    src/solve.f)