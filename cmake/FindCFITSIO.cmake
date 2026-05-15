# FindCFITSIO.cmake - locate CFITSIO library
# Sets: CFITSIO_FOUND, CFITSIO_INCLUDE_DIRS, CFITSIO_LIBRARIES
# Creates imported target: CFITSIO::CFITSIO

find_path(CFITSIO_INCLUDE_DIR
    NAMES fitsio.h
    PATHS
        $ENV{CONDA_PREFIX}/Library/include
        $ENV{CONDA_PREFIX}/include
        /usr/include
        /usr/local/include
    PATH_SUFFIXES cfitsio
)

find_library(CFITSIO_LIBRARY
    NAMES cfitsio libcfitsio
    PATHS
        $ENV{CONDA_PREFIX}/Library/lib
        $ENV{CONDA_PREFIX}/lib
        /usr/lib
        /usr/local/lib
        /usr/lib/x86_64-linux-gnu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFITSIO
    REQUIRED_VARS CFITSIO_LIBRARY CFITSIO_INCLUDE_DIR
)

if(CFITSIO_FOUND AND NOT TARGET CFITSIO::CFITSIO)
    add_library(CFITSIO::CFITSIO UNKNOWN IMPORTED)
    set_target_properties(CFITSIO::CFITSIO PROPERTIES
        IMPORTED_LOCATION "${CFITSIO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CFITSIO_INCLUDE_DIR}"
    )
endif()

set(CFITSIO_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY})
mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)
