# phosimcatgen executable - star/galaxy catalog generator
#
# The original Makefile hard-codes CC = g++; CMake uses
# ${CMAKE_CXX_COMPILER} automatically.

add_executable(phosimcatgen
    phosimcatgen.cpp
    cgparams.cpp
    catalogobject.cpp
    starobject.cpp
    galaxyobject.cpp
    galaxycomplex.cpp
    sphericalcell.cpp
)

target_include_directories(phosimcatgen PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(phosimcatgen PRIVATE
    phosim_ancillary
    CFITSIO::CFITSIO
)
