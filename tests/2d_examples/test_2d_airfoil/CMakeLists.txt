CMAKE_MINIMUM_REQUIRED(VERSION 3.8.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake) # main (top) cmake dir

set(CMAKE_VERBOSE_MAKEFILE on)

STRING( REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR} )
PROJECT("${CURRENT_FOLDER}")
include(ImportSPHINXsysFromSource_for_2D_build)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/airfoil_flap_front.dat
        DESTINATION ${BUILD_INPUT_PATH})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/airfoil_flap_rear.dat
        DESTINATION ${BUILD_INPUT_PATH})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/airfoil_wing.dat
        DESTINATION ${BUILD_INPUT_PATH})

aux_source_directory(. DIR_SRCS)
ADD_EXECUTABLE(${PROJECT_NAME} ${EXECUTABLE_OUTPUT_PATH} ${DIR_SRCS} )

add_test(NAME ${PROJECT_NAME}_particle_relaxation 
		 COMMAND ${PROJECT_NAME} --r=true
                 WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
                 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
    target_link_libraries(${PROJECT_NAME} sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES})
    add_dependencies(${PROJECT_NAME} sphinxsys_2d sphinxsys_static_2d)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    	target_link_libraries(${PROJECT_NAME} sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} ${Boost_LIBRARIES} stdc++)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(${PROJECT_NAME} sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES}  ${Boost_LIBRARIES} stdc++ stdc++fs)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
