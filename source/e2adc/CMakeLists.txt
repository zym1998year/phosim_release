# e2adc executable - electron-to-ADU converter

add_executable(e2adc
    e2adc.cpp
)

target_include_directories(e2adc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(e2adc PRIVATE
    phosim_ancillary
    CFITSIO::CFITSIO
)
