# Generate bin/version file mimicking original Makefile behavior

set(PHOSIM_VERSION_FILE ${CMAKE_SOURCE_DIR}/bin/version)

find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} log -1 --format=%H
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE PHOSIM_GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE PHOSIM_GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT PHOSIM_GIT_COMMIT)
        set(PHOSIM_GIT_COMMIT "none")
    endif()
    if(NOT PHOSIM_GIT_DESCRIBE)
        set(PHOSIM_GIT_DESCRIBE "unknown")
    endif()
else()
    set(PHOSIM_GIT_COMMIT "none")
    set(PHOSIM_GIT_DESCRIBE "${CMAKE_SOURCE_DIR}")
endif()

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(WRITE ${PHOSIM_VERSION_FILE} "commit ${PHOSIM_GIT_COMMIT}\n${PHOSIM_GIT_DESCRIBE}\n")
