# trim executable - PhoSim's input catalog trimming utility

add_executable(trim
    main.cpp
    trim.cpp
)

target_include_directories(trim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/source
)

target_link_libraries(trim PRIVATE
    phosim_ancillary
    ZLIB::ZLIB
)

# Test executable
# test_trim.cpp uses #include "validation/unittest.h" (relative to repo root)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_trim.cpp)
    add_executable(test_trim
        test_trim.cpp
        trim.cpp
    )
    target_include_directories(test_trim PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/source
        ${CMAKE_SOURCE_DIR}
    )
    target_link_libraries(test_trim PRIVATE phosim_ancillary ZLIB::ZLIB)
    add_test(NAME trim_smoke COMMAND test_trim)
endif()
