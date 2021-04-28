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

#FOREACH(file1 ${SCR_FILES})
	#message(STATUS ${file1})
#ENDFOREACH()

ADD_LIBRARY(sphinxsys_3d SHARED ${SCR_FILES})
ADD_LIBRARY(sphinxsys_static_3d STATIC ${SCR_FILES})

SET_TARGET_PROPERTIES(sphinxsys_static_3d PROPERTIES OUTPUT_NAME "sphinxsys_3d")
#SET_TARGET_PROPERTIES(sphinxsys PROPERTIES VERSION 1.0 SOVERSION 0)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    target_link_libraries(sphinxsys_3d ${TBB_LIBRARYS} ${Simbody_LIBRARIES})
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    	target_link_libraries(sphinxsys_3d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} ${Boost_LIBRARIES} stdc++)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(sphinxsys_3d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} stdc++ stdc++fs)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

INSTALL(TARGETS sphinxsys_3d sphinxsys_static_3d
RUNTIME DESTINATION 3d_code/bin
LIBRARY DESTINATION 3d_code/lib
ARCHIVE DESTINATION 3d_code/lib)

FILE(GLOB_RECURSE hpp_headers  ${PROJECT_SOURCE_DIR}/src/shared/*.hpp ${PROJECT_SOURCE_DIR}/src/for_3D_build/*.hpp)
INSTALL(FILES ${hpp_headers} DESTINATION 3d_code/include)
