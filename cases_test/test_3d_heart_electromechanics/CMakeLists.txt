CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake) # main (top) cmake dir
set(CMAKE_VERBOSE_MAKEFILE on)

STRING( REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR} )
PROJECT("${CURRENT_FOLDER}")

include(ImportSPHINXsysFromSource_for_3D_build)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

if(MSVC)
	SET(BUILD_INPUT_PATH "${PROJECT_BINARY_DIR}/input")
	SET(BUILD_RELOAD_PATH "${PROJECT_BINARY_DIR}/reload")
else(MSVC)
	SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
	SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")
endif(MSVC)

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_INPUT_PATH})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/heart-new.stl
        DESTINATION ${BUILD_INPUT_PATH})

aux_source_directory(. DIR_SRCS)
ADD_EXECUTABLE(${PROJECT_NAME} ${DIR_SRCS})

add_test(NAME ${PROJECT_NAME}_particle_relaxation 
		 COMMAND ${PROJECT_NAME} --r=true
		 WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
add_test(NAME ${PROJECT_NAME} 
		 COMMAND ${PROJECT_NAME}  --i=true
		 WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(${PROJECT_NAME} sphinxsys_3d ${TBB_LIBRARYS} debug ${Simbody_DEBUG_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} sphinxsys_3d ${TBB_LIBRARYS} optimized ${Simbody_RELEASE_LIBRARIES})
    add_dependencies(${PROJECT_NAME} sphinxsys_3d sphinxsys_static_3d)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    	target_link_libraries(${PROJECT_NAME} sphinxsys_3d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} ${Boost_LIBRARIES} stdc++)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(${PROJECT_NAME} sphinxsys_3d ${TBB_LIBRARYS} ${Simbody_LIBRARIES}  ${Boost_LIBRARIES} stdc++ stdc++fs)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")