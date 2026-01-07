cmake_minimum_required(VERSION 3.12)

project(bode
    VERSION 0.1.0
    LANGUAGES Fortran C)

add_library(bode
	src/bode_module.F90
	src/bode_nonlin.f90
	src/bode_solve.f90)

add_executable(test_banded1 test/test_banded1.f90)
target_link_libraries(test_banded1 bode)

add_executable(test_banded2 test/test_banded2.f90)
target_link_libraries(test_banded2 bode)

add_executable(robertson examples/robertson.F90)
target_link_libraries(robertson bode)
target_link_libraries(robertson lapack)

add_executable(robertson_old examples/robertson.F90 src/bode.f)
target_compile_definitions(robertson_old PRIVATE ORIG)

add_executable(test_vdp examples/test_van_der_pol.f90)
target_link_libraries(test_vdp bode)

add_executable(test_vdp_old examples/test_van_der_pol.f90 src/bode.f)

# Banded procedures by Jeff Thorson
add_library(bandfs bandfs/bandf.f bandfs/bands.f bandfs/bndmv.f90 bandfs/bndsl.f90)

