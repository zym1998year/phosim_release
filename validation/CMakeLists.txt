# phosim_validation executable - PhoSim internal unit tests
#
# unittest.cpp uses functions from raytrace.cpp (normalize, reflect,
# refract, propagate) and inline helpers from helpers.h.  raytrace.cpp
# has no main(); its entry point is in raytrace/main.cpp, so there is
# no symbol collision.  We add raytrace.cpp directly to SOURCES rather
# than re-linking the full raytrace executable.
#
# The original Makefile only listed fftw3 as an external lib; raytrace.cpp
# itself only needs math (already pulled via phosim_ancillary PUBLIC links).

add_executable(phosim_validation
    unittest.cpp
    ${CMAKE_SOURCE_DIR}/source/raytrace/raytrace.cpp
)

target_include_directories(phosim_validation PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
    ${CMAKE_SOURCE_DIR}/source/raytrace
)

target_link_libraries(phosim_validation PRIVATE
    phosim_ancillary
    FFTW3::fftw3
    Threads::Threads
)

add_test(NAME validation_unittest COMMAND phosim_validation)
