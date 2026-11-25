STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
# Reuse CPU case STL / CSV / JSON (already meters; do not rescale).
set(SPDB_PROG_TEMP_INPUT
    ${CMAKE_SOURCE_DIR}/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/input)
file(COPY ${SPDB_PROG_TEMP_INPUT}/ DESTINATION ${BUILD_INPUT_PATH})

aux_source_directory(. DIR_SRCS)
add_executable(${PROJECT_NAME} ${DIR_SRCS})
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
target_link_libraries(${PROJECT_NAME} extra_sources_3d)

add_test(NAME ${PROJECT_NAME}_geometry COMMAND ${PROJECT_NAME} --geometry_only
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
set_tests_properties(${PROJECT_NAME}_geometry PROPERTIES LABELS "spdb;sycl")
