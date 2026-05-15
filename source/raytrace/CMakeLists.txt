# raytrace executable - PhoSim's main ray tracing engine
#
# NOTE on #include-d cpp files:
# atmospheresetup.cpp, telescopesetup.cpp, photonmanipulate.cpp,
# photonoptimization.cpp, photonloop.cpp, particleloop.cpp, sourceloop.cpp,
# coordinate.cpp, settings.cpp, sed.cpp, ext.cpp are all pulled in via
# #include inside image.cpp / observation.cpp and must NOT be listed as
# separate compilation units here.
#
# NOTE on include paths:
# main.cpp uses #include "raytrace/image.h", and raytrace sources use
# #include "ancillary/random.h" - both resolve relative to source/.
# We therefore add CMAKE_SOURCE_DIR/source as a private include root.

set(RAYTRACE_SOURCES
    main.cpp
    observation.cpp
    image.cpp
    raytrace.cpp
    galaxy.cpp
    dust.cpp
    event.cpp
    grating.cpp
    counter.cpp
    silicon.cpp
    surface.cpp
    coating.cpp
    perturbation.cpp
    obstruction.cpp
    screen.cpp
    medium.cpp
    air.cpp
    chip.cpp
    contamination.cpp
    fea.cpp
    mie.cpp
)

add_executable(raytrace ${RAYTRACE_SOURCES})

target_include_directories(raytrace PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(raytrace PRIVATE
    phosim_ancillary
    CFITSIO::CFITSIO
    FFTW3::fftw3
    ZLIB::ZLIB
    Threads::Threads
)
