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
# Copy COMSOL TEAM7 reference tables (if present) for curve validation.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/input/comsol_team7_raw)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/comsol_team7_raw/ DESTINATION ${BUILD_INPUT_PATH}/comsol_team7_raw)
endif()

add_executable(${PROJECT_NAME} particle_generation_em.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE sphinxsys_3d extra_sources_3d)
# Run with working directory = bin/ so ./input/coil.stl and ./input/plate.stl are found
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
