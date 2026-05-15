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

target_link_libraries(phosim_ancillary PUBLIC
    CFITSIO::CFITSIO
    FFTW3::fftw3
    ZLIB::ZLIB
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
endif()
