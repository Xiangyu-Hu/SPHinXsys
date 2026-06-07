# Sync TEAM7 native SI reload + STL input into test bin/ at configure time.
# Priority:
#   1) TEAM7_NATIVE_RELOAD_DIR (CMake cache)
#   2) ${CMAKE_BINARY_DIR}/.../particle_generation_em/bin/reload (after relax in same build tree)
#   3) ${CMAKE_CURRENT_SOURCE_DIR}/data/reload (optional committed snapshot)

if(NOT DEFINED BUILD_RELOAD_PATH OR NOT DEFINED BUILD_INPUT_PATH)
    message(FATAL_ERROR "team7_native_reload_sync.cmake requires BUILD_RELOAD_PATH and BUILD_INPUT_PATH")
endif()

if(NOT DEFINED TEAM7_NATIVE_RELOAD_DIR)
    set(TEAM7_NATIVE_RELOAD_DIR "" CACHE PATH "Directory with Coil_rld.xml, Plate_rld.xml, Air_rld.xml (SI TEAM7)")
endif()

set(_TEAM7_RELOAD_SRC "")
if(TEAM7_NATIVE_RELOAD_DIR AND EXISTS "${TEAM7_NATIVE_RELOAD_DIR}")
    set(_TEAM7_RELOAD_SRC "${TEAM7_NATIVE_RELOAD_DIR}")
elseif(EXISTS "${CMAKE_BINARY_DIR}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin/reload")
    set(_TEAM7_RELOAD_SRC
        "${CMAKE_BINARY_DIR}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin/reload")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/data/reload")
    set(_TEAM7_RELOAD_SRC "${CMAKE_CURRENT_SOURCE_DIR}/data/reload")
endif()

if(_TEAM7_RELOAD_SRC)
    file(COPY "${_TEAM7_RELOAD_SRC}/" DESTINATION "${BUILD_RELOAD_PATH}")
    message(STATUS "TEAM7 native reload: copied from ${_TEAM7_RELOAD_SRC}")
else()
    message(WARNING "TEAM7 native reload: no reload XML found. Run particle_generation_em or set -DTEAM7_NATIVE_RELOAD_DIR=...")
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../particle_generation_em/data)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../particle_generation_em/data/ DESTINATION ${BUILD_INPUT_PATH})
elseif(EXISTS "${CMAKE_BINARY_DIR}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin/input")
    file(COPY "${CMAKE_BINARY_DIR}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin/input/"
         DESTINATION ${BUILD_INPUT_PATH})
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../reference_data/team7)
    if(NOT DEFINED BUILD_REFERENCE_PATH)
        set(BUILD_REFERENCE_PATH "${EXECUTABLE_OUTPUT_PATH}/reference_data/team7")
        file(MAKE_DIRECTORY ${BUILD_REFERENCE_PATH})
    endif()
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../reference_data/team7/ DESTINATION ${BUILD_REFERENCE_PATH})
endif()
