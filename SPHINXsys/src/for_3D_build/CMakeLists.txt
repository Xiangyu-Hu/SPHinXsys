## 3D build
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake) # main (top) cmake dir
include(Dirsearch_for_3D_build)

## prepare dirctoriesfor head and source files
HEADER_DIRECTORIES_SHARED(headdirs_shared)
SOURCE_DIRECTORIES_SHARED(sourcedirs_shared)
HEADER_DIRECTORIES_3D(headdirs_3D)
SOURCE_DIRECTORIES_3D(sourcedirs_3D)

SET(usefuldirs  ${headdirs_3D} ${sourcedirs_3D})
LIST(REMOVE_DUPLICATES usefuldirs)

SET(usefulsubdirs ${usefuldirs})
LIST(REMOVE_ITEM usefulsubdirs ${CMAKE_CURRENT_SOURCE_DIR})

#Add all useful subdirectories
FOREACH(subdir_path ${usefulsubdirs})
	#message(STATUS ${subdir_path})
	ADD_SUBDIRECTORY(${subdir_path})
ENDFOREACH()

# combin head and souce directories
SET(headdirs ${headdirs_shared} ${headdirs_3D})
SET(sourcedirs ${sourcedirs_shared} ${sourcedirs_3D})

##Add all header dirs
FOREACH(headdir_path ${headdirs})
	#message(STATUS ${headdir_path})
	INCLUDE_DIRECTORIES("${headdir_path}")
ENDFOREACH()

##Add all source files
set(SCR_FILES "")
FOREACH(srcdir_path ${sourcedirs})
	#message(STATUS ${srcdir_path})
	set(DIR_scrs "")
	AUX_SOURCE_DIRECTORY(${srcdir_path} DIR_scrs)
	list(APPEND SCR_FILES ${DIR_scrs})
ENDFOREACH()

if(NOT SPH_ONLY_STATIC_BUILD)
	### SPHinXsys dynamic lib ###
	ADD_LIBRARY(sphinxsys_3d SHARED ${SCR_FILES})

	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		target_link_libraries(sphinxsys_3d)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(sphinxsys_3d stdc++)
		else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			target_link_libraries(sphinxsys_3d stdc++ stdc++fs)
		endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

		if(DEFINED BOOST_AVAILABLE) # link Boost if available (not for Windows)
			target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
		endif()
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

	if(NOT BUILD_WITH_ONETBB) # link TBB if not built by the project
		target_link_libraries(sphinxsys_3d ${TBB_LIBRARYS}})
	else()
		target_link_libraries(sphinxsys_3d TBB::tbb TBB::tbbmalloc) # TBB::tbbmalloc_proxy is not needed
	endif()
	if(NOT BUILD_WITH_SIMBODY) # link Simbody if not built by the project
		target_link_libraries(sphinxsys_3d ${Simbody_LIBRARIES})
	endif()
	### SPHinXsys dynamic lib ###
endif()

### SPHinXsys static lib ###
ADD_LIBRARY(sphinxsys_static_3d STATIC ${SCR_FILES})

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	target_link_libraries(sphinxsys_3d)
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(sphinxsys_3d stdc++)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(sphinxsys_3d stdc++ stdc++fs)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

	if(DEFINED BOOST_AVAILABLE) # link Boost if available (not for Windows)
		target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
	endif()
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(NOT BUILD_WITH_ONETBB) # link TBB if not built by the project
	target_link_libraries(sphinxsys_static_3d ${TBB_LIBRARYS}})
else()
	target_link_libraries(sphinxsys_static_3d TBB::tbb TBB::tbbmalloc) # TBB::tbbmalloc_proxy is not needed
endif()
if(NOT BUILD_WITH_SIMBODY) # link Simbody if not built by the project
	target_link_libraries(sphinxsys_static_3d ${Simbody_LIBRARIES})
endif()
### SPHinXsys static lib ###

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

if(NOT SPH_ONLY_STATIC_BUILD)
	INSTALL(TARGETS sphinxsys_3d sphinxsys_static_3d
	RUNTIME DESTINATION 3d_code/bin
	LIBRARY DESTINATION 3d_code/lib
	ARCHIVE DESTINATION 3d_code/lib)
endif()

INSTALL(TARGETS sphinxsys_static_3d
RUNTIME DESTINATION 3d_code/bin
LIBRARY DESTINATION 3d_code/lib
ARCHIVE DESTINATION 3d_code/lib)

# FILE(GLOB_RECURSE hpp_headers  ${PROJECT_SOURCE_DIR}/src/shared/*.hpp ${PROJECT_SOURCE_DIR}/src/for_3D_build/*.hpp)
# INSTALL(FILES ${hpp_headers} DESTINATION 3d_code/include)
