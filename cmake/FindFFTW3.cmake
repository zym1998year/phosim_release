# FindFFTW3.cmake - locate FFTW3 library
# Sets: FFTW3_FOUND, FFTW3_INCLUDE_DIRS, FFTW3_LIBRARIES
# Creates imported target: FFTW3::fftw3

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    PATHS
        $ENV{CONDA_PREFIX}/Library/include
        $ENV{CONDA_PREFIX}/include
        /usr/include
        /usr/local/include
)

find_library(FFTW3_LIBRARY
    NAMES fftw3 libfftw3 libfftw3-3
    PATHS
        $ENV{CONDA_PREFIX}/Library/lib
        $ENV{CONDA_PREFIX}/lib
        /usr/lib
        /usr/local/lib
        /usr/lib/x86_64-linux-gnu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
)

if(FFTW3_FOUND AND NOT TARGET FFTW3::fftw3)
    add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
    set_target_properties(FFTW3::fftw3 PROPERTIES
        IMPORTED_LOCATION "${FFTW3_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
    )
endif()

set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)
