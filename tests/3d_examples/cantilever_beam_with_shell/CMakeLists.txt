set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SPHINXSYS_PROJECT_DIR}/cmake) # main (top) cmake dir
set(CMAKE_VERBOSE_MAKEFILE on)

STRING( REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR} )
PROJECT("${CURRENT_FOLDER}")

include(ImportSPHINXsysFromSource_for_3D_build)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
set(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
set(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

aux_source_directory(. DIR_SRCS)
ADD_EXECUTABLE(${PROJECT_NAME} ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME} gtest gtest_main)

#gtest_discover_tests(${PROJECT_NAME} WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

if(NOT SPH_ONLY_STATIC_BUILD) # usual dynamic build
	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
		target_link_libraries(${PROJECT_NAME} sphinxsys_3d)
		add_dependencies(${PROJECT_NAME} sphinxsys_3d)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(${PROJECT_NAME} sphinxsys_3d stdc++)
		else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(${PROJECT_NAME} sphinxsys_3d stdc++ stdc++fs dl)
		endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

		if(DEFINED BOOST_AVAILABLE) # link Boost if available (not for Windows)
			target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
		endif()
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
else() # static build only
	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
                                set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")
		target_link_libraries(${PROJECT_NAME} sphinxsys_static_3d)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(${PROJECT_NAME} sphinxsys_static_3d stdc++)
		else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(${PROJECT_NAME} sphinxsys_static_3d stdc++ stdc++fs dl)
		endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

		if(DEFINED BOOST_AVAILABLE) # link Boost if available (not for Windows)
			target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
		endif()
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
endif()
if(NOT BUILD_WITH_SIMBODY) # link Simbody if not built by the project
target_link_libraries(${PROJECT_NAME} ${Simbody_LIBRARIES})
endif()
if(NOT BUILD_WITH_ONETBB) # link TBB if not built by the project
target_link_libraries(${PROJECT_NAME} ${TBB_LIBRARYS})
endif()