set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake) # main (top) cmake dir

set(CMAKE_VERBOSE_MAKEFILE on)

include(ImportSPHINXsysFromSource_for_2D_build)

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

aux_source_directory(. DIR_SRCS)
ADD_EXECUTABLE(${PROJECT_NAME} ${DIR_SRCS})

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(${PROJECT_NAME} sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES})
    add_dependencies(${PROJECT_NAME} sphinxsys_2d)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    	target_link_libraries(${PROJECT_NAME} sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} ${Boost_LIBRARIES} stdc++)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(${PROJECT_NAME} sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES}  ${Boost_LIBRARIES} stdc++ stdc++fs)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
