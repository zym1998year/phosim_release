# instrument executable - telescope/camera modeling

add_executable(instrument
    instrument.cpp
    instrumentfiles.cpp
    distortion.cpp
)

target_include_directories(instrument PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(instrument PRIVATE
    phosim_ancillary
    Threads::Threads
)
