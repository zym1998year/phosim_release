cmake_minimum_required(VERSION 3.20)
project(phosim VERSION 6.3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(PHOSIM_BUILD_VENDORED_DEPS "Build cfitsio/fftw from vendored tarballs (Linux only)" OFF)
option(PHOSIM_EXT "Enable EXT mode" OFF)
option(PHOSIM_MSVC_FAST_FLOAT "Use /fp:fast on MSVC (default OFF; off main CI)" OFF)

if(WIN32 AND PHOSIM_BUILD_VENDORED_DEPS)
    message(FATAL_ERROR "PHOSIM_BUILD_VENDORED_DEPS=ON is not supported on Windows. Use conda-forge cfitsio/fftw instead.")
endif()

# Compiler flags
if(MSVC)
    add_compile_options(/permissive- /EHsc /Zc:__cplusplus /utf-8 /DNOMINMAX)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS NOMINMAX _USE_MATH_DEFINES)
    # PhoSim uses large fixed stack arrays (e.g. double[10000] in instrument.cpp) that
    # fit Linux's 8 MB default stack but overflow MSVC's 1 MB default. Match Linux (8 MB).
    add_link_options(/STACK:8388608)
    if(PHOSIM_MSVC_FAST_FLOAT)
        add_compile_options($<$<CONFIG:Release>:/fp:fast>)
    else()
        add_compile_options($<$<CONFIG:Release>:/fp:precise>)
    endif()
    add_compile_options($<$<CONFIG:Release>:/O2>)
else()
    add_compile_options(-Wall -fvisibility=default)
    # Note: each flag must be a separate argument; do NOT bundle multiple flags
    # into a single generator expression (they would be passed as one token).
    add_compile_options(
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-ffast-math>
        $<$<CONFIG:Release>:-fno-finite-math-only>
    )
endif()

if(PHOSIM_BUILD_VENDORED_DEPS)
    message(FATAL_ERROR "Vendored deps build is a follow-up; use conda-forge or system packages for now.")
endif()

find_package(CFITSIO REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Output directories (single-config + multi-config compatible)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
foreach(CONFIG IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER ${CONFIG} CONFIG_U)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CONFIG_U} ${CMAKE_SOURCE_DIR}/bin)
endforeach()

# Version file generation
include(PhosimVersion)

# ext.h generation (replaces grep on bin/setup)
configure_file(
    ${CMAKE_SOURCE_DIR}/source/raytrace/ext.h.in
    ${CMAKE_BINARY_DIR}/generated/source/raytrace/ext.h
    @ONLY
)
include_directories(${CMAKE_BINARY_DIR}/generated/source/raytrace)

# phosim_compat.h generation (PHOSIM_NOINLINE and other compiler-portability macros)
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/phosim_compat.h.in
    ${CMAKE_BINARY_DIR}/generated/phosim_compat.h
    @ONLY
)
include_directories(${CMAKE_BINARY_DIR}/generated)

# CTest must be enabled before any add_subdirectory that calls add_test()
enable_testing()

# Subdirectories - all modules enabled
add_subdirectory(source/ancillary)
add_subdirectory(source/raytrace)
add_subdirectory(source/trim)
add_subdirectory(source/e2adc)
add_subdirectory(source/instrument)
add_subdirectory(source/atmosphere)
add_subdirectory(source/catalogtools)
add_subdirectory(validation)
