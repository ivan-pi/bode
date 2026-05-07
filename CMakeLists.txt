cmake_minimum_required(VERSION 3.18)

project(bode VERSION 0.2.0 LANGUAGES Fortran C)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

option(BODE_ENABLE_WARNINGS "Enable compiler warnings" ON)
option(BODE_BUILD_EXAMPLES "Build example programs" ON)
option(BODE_BUILD_LEGACY "Build legacy comparison executables" OFF)
option(BODE_BUILD_DOCS "Enable FORD documentation target" OFF)

find_package(LAPACK REQUIRED)

add_library(bode
    src/bode_module.F90
    src/bode_nonlin.f90
    src/bode_solve.f90
)
add_library(bode::bode ALIAS bode)

target_link_libraries(bode PUBLIC LAPACK::LAPACK)

target_include_directories(bode
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(bandfs
    bandfs/bandf.f
    bandfs/bands.f
    bandfs/bndmv.f90
    bandfs/bndsl.f90
)
add_library(bode::bandfs ALIAS bandfs)

target_include_directories(bandfs
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bandfs/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(BODE_ENABLE_WARNINGS)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Intel")
        target_compile_options(bode PRIVATE -Wall -Wextra -Wimplicit-interface -fcheck=all)
        target_compile_options(bandfs PRIVATE -Wall -Wextra)
    endif()
endif()

if(BUILD_TESTING)
    add_library(test_support test/test_support.f90)
    target_link_libraries(test_support PUBLIC bode)

    add_executable(test_banded1 test/test_banded1.f90)
    target_link_libraries(test_banded1 PRIVATE test_support bode)

    add_executable(test_banded2 test/test_banded2.f90)
    target_link_libraries(test_banded2 PRIVATE test_support bode bandfs)

    add_executable(test_bandfs_api test/test_bandfs_api.f90)
    target_link_libraries(test_bandfs_api PRIVATE test_support bode bandfs)

    add_executable(test_vdp_regression test/test_vdp_regression.f90)
    target_link_libraries(test_vdp_regression PRIVATE test_support bode)

    add_executable(test_robertson_regression test/test_robertson_regression.f90)
    target_link_libraries(test_robertson_regression PRIVATE test_support bode)

    add_executable(test_invalid_iopt test/test_invalid_iopt.f90)
    target_link_libraries(test_invalid_iopt PRIVATE bode)

    add_test(NAME unit.banded_ltri_tsol COMMAND test_banded1)
    set_tests_properties(unit.banded_ltri_tsol PROPERTIES LABELS "unit")

    add_test(NAME unit.bandfs_api COMMAND test_bandfs_api)
    set_tests_properties(unit.bandfs_api PROPERTIES LABELS "unit")

    add_test(NAME regression.banded2 COMMAND test_banded2
        ${CMAKE_CURRENT_SOURCE_DIR}/test/data/banded2_reference.csv)
    set_tests_properties(regression.banded2 PROPERTIES LABELS "regression")

    add_test(NAME regression.vdp COMMAND test_vdp_regression
        ${CMAKE_CURRENT_SOURCE_DIR}/test/data/vdp_reference.txt)
    set_tests_properties(regression.vdp PROPERTIES LABELS "regression")

    add_test(NAME regression.robertson COMMAND test_robertson_regression
        ${CMAKE_CURRENT_SOURCE_DIR}/test/data/robertson_reference.txt)
    set_tests_properties(regression.robertson PROPERTIES LABELS "regression")

    add_test(NAME failure.invalid_iopt COMMAND test_invalid_iopt)
    set_tests_properties(failure.invalid_iopt PROPERTIES LABELS "failure" WILL_FAIL TRUE)
endif()

if(BODE_BUILD_EXAMPLES)
    add_executable(example_robertson examples/robertson.F90)
    target_link_libraries(example_robertson PRIVATE bode)

    add_executable(example_vdp examples/test_van_der_pol.f90)
    target_link_libraries(example_vdp PRIVATE bode)

    add_executable(example_getting_started examples/getting_started.f90)
    target_link_libraries(example_getting_started PRIVATE bode)

    add_executable(example_dense_vs_banded examples/dense_vs_banded_robertson.f90)
    target_link_libraries(example_dense_vs_banded PRIVATE bode)

    if(BODE_BUILD_LEGACY)
        add_executable(robertson_old examples/robertson.F90 src/bode.f)
        target_compile_definitions(robertson_old PRIVATE ORIG)

        add_executable(test_vdp_old examples/test_van_der_pol.f90 src/bode.f)
    endif()
endif()

if(BODE_BUILD_DOCS)
    find_program(FORD_EXECUTABLE ford)
    if(FORD_EXECUTABLE)
        add_custom_target(docs
            COMMAND ${FORD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ford.md
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating FORD API documentation"
            VERBATIM
        )
    else()
        message(WARNING "FORD not found; docs target will be unavailable")
    endif()
endif()

install(TARGETS bode bandfs
    EXPORT bodeTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES
    bandfs/include/bandfs.fi
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT bodeTargets
    FILE bodeTargets.cmake
    NAMESPACE bode::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bode
)

configure_package_config_file(
    cmake/bodeConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/bodeConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bode
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/bodeConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/bodeConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/bodeConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bode
)

set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/\\.git/"
    "\\.swp$"
    "\\.DS_Store"
)
include(CPack)
