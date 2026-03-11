STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
# Copy data/ (coil.stl, plate.stl) to input/ so ./input/coil.stl is found when running from bin/
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/data)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION ${BUILD_INPUT_PATH})
endif()

add_executable(${PROJECT_NAME})
aux_source_directory(. DIR_SRCS)
target_sources(${PROJECT_NAME} PRIVATE ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME} sphinxsys_3d)
# Run with working directory = bin/ so ./input/coil.stl and ./input/plate.stl are found
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
