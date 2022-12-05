set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake) # main (top) cmake dir
set(CMAKE_VERBOSE_MAKEFILE on)

STRING( REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR} )
PROJECT("${CURRENT_FOLDER}")

include(ImportSPHINXsysFromSource_for_3D_build)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
set(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
set(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})

set(FILES_STL	"bernoulli_beam_20x.stl")
foreach(STL_FILE ${FILES_STL})
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/${STL_FILE} DESTINATION ${BUILD_INPUT_PATH})
endforeach()

aux_source_directory(. DIR_SRCS)
ADD_EXECUTABLE(${PROJECT_NAME} ${EXECUTABLE_OUTPUT_PATH} ${DIR_SRCS})
target_link_directories(${PROJECT_NAME} PRIVATE ${TBBLIBDIR})
target_link_libraries(${PROJECT_NAME} structural_simulation_module ${TBBLIB})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/SPHinXsys/modules/structural_simulation ${SPHINXSYS_PROJECT_DIR}/3rd_party/eigen)

if(NOT EMSCRIPTEN)
	add_test(NAME ${PROJECT_NAME} 
			COMMAND ${PROJECT_NAME} --r=true
			WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
endif()