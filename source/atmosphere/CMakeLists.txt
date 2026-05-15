# atmosphere executable - atmospheric turbulence/seeing model
#
# NOTE: atmosphere.cpp #includes atmosphere_creator.cpp, turb2d.cpp,
# cloud.cpp, airglow.cpp internally; main.cpp #includes operator.cpp.
# These files must NOT be listed as separate compilation units.
# air.cpp from raytrace/ is a separate object (no internal #include),
# so it must be listed explicitly here, mirroring the original Makefile.

add_executable(atmosphere
    main.cpp
    atmosphere.cpp
    ${CMAKE_SOURCE_DIR}/source/raytrace/air.cpp
)

target_include_directories(atmosphere PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
    ${CMAKE_SOURCE_DIR}/source/raytrace
)

target_link_libraries(atmosphere PRIVATE
    phosim_ancillary
    CFITSIO::CFITSIO
    FFTW3::fftw3
    Threads::Threads
)

# test_atmosphere.cpp uses non-standard 'void main()' which GCC rejects
# as an error. This test is excluded from the build until the source is
# fixed in a later task (C++ source corrections are out of scope for Task 1).
# if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_atmosphere.cpp)
#     add_executable(test_atmosphere test_atmosphere.cpp)
#     ...
# endif()
