# phosim_ancillary - shared sources used by raytrace, trim, e2adc, etc.

add_library(phosim_ancillary STATIC
    random.cpp
    readtext.cpp
)

target_include_directories(phosim_ancillary PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    # readtext.cpp uses #include "ancillary/readtext.h" which resolves
    # relative to source/; expose source/ as a public include root so that
    # consumers (raytrace etc.) can also use the same path.
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(phosim_ancillary
    # fits.h exposes fitsfile* in every function signature -> PUBLIC
    PUBLIC CFITSIO::CFITSIO
    # fft.h exposes fftw_complex / fftw_plan in function signatures -> PUBLIC
    PUBLIC FFTW3::fftw3
    # readtext.h does not expose any zlib types; gzip usage is an
    # implementation detail in readtext.cpp only -> PRIVATE
    PRIVATE ZLIB::ZLIB
)

# Tests
# test_random.cpp and test_readtext.cpp both use #include "validation/unittest.h"
# which resolves relative to the repo root; add CMAKE_SOURCE_DIR to include path.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_random.cpp)
    add_executable(test_random test_random.cpp)
    target_include_directories(test_random PRIVATE ${CMAKE_SOURCE_DIR})
    target_link_libraries(test_random PRIVATE phosim_ancillary)
    add_test(NAME ancillary_random COMMAND test_random)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_readtext.cpp)
    add_executable(test_readtext test_readtext.cpp)
    target_include_directories(test_readtext PRIVATE ${CMAKE_SOURCE_DIR})
    target_link_libraries(test_readtext PRIVATE phosim_ancillary)
    add_test(NAME ancillary_readtext COMMAND test_readtext)
    # NOTE: test_readtext reads parameter data from std::cin; CTest cannot
    # pipe stdin without a shell wrapper, which is non-portable on Windows.
    # The binary itself builds and links correctly (test_random exercises the
    # Random class; ancillary_readtext tests the readText class at the
    # filesystem level).  The fixture files exist in source/ancillary/
    # (test.txt) and data/atmosphere/ (nprofile.txt, o2cs.txt).  A proper
    # integration test will be wired once the POSIX/portability layer is
    # settled in Task 3+.  For now mark disabled so ctest exits cleanly.
    set_tests_properties(ancillary_readtext PROPERTIES
        DISABLED TRUE
    )
endif()
