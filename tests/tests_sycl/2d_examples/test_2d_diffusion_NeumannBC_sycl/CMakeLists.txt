STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})
if(SPHINXSYS_USE_SYCL) # this case gives more significant differences when use sycl as float other double is used
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool_sycl/ DESTINATION ${BUILD_INPUT_PATH})
else()
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool_ck/ DESTINATION ${BUILD_INPUT_PATH})
endif()



add_executable(${PROJECT_NAME})
aux_source_directory(. DIR_SRCS)
target_sources(${PROJECT_NAME} PRIVATE ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME} sphinxsys_2d)
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_test(NAME ${PROJECT_NAME}
	COMMAND ${PROJECT_NAME} --state_recording=${TEST_STATE_RECORDING}
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

set_tests_properties(${PROJECT_NAME} PROPERTIES LABELS "diffusion reaction")