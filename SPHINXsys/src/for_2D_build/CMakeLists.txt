## 2D build
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake) # main (top) cmake dir

## prepare dirctories for head and source files
HEADER_DIRECTORIES_SHARED(headdirs_shared)
SOURCE_DIRECTORIES_SHARED(sourcedirs_shared)

include(Dirsearch_for_2D_build)
HEADER_DIRECTORIES_2D(headdirs_2D)
SOURCE_DIRECTORIES_2D(sourcedirs_2D)

SET(usefuldirs ${headdirs_2D} ${sourcedirs_2D})
LIST(REMOVE_DUPLICATES usefuldirs)

SET(usefulsubdirs ${usefuldirs})
LIST(REMOVE_ITEM usefulsubdirs ${CMAKE_CURRENT_SOURCE_DIR})

##Add all useful subdirectories
FOREACH(subdir_path ${usefulsubdirs})
	#message(STATUS ${subdir_path})
	ADD_SUBDIRECTORY(${subdir_path})
ENDFOREACH()

if(BUILD_WITH_IMAGE_PROCESS)
	include(Dirsearch_for_image_process)
	HEADER_DIRECTORIES_IMAGE_PROCESS(headdirs_image_process)
	SOURCE_DIRECTORIES_IMAGE_PROCESS(sourcedirs_image_process)
	## combin head and souce directories
	SET(headdirs ${headdirs_shared} ${headdirs_2D}  ${headdirs_image_process})
	SET(sourcedirs ${sourcedirs_shared} ${sourcedirs_2D} ${sourcedirs_image_process})
else(BUILD_WITH_IMAGE_PROCESS)
	## combin head and souce directories
	SET(headdirs ${headdirs_shared} ${headdirs_2D})
	SET(sourcedirs ${sourcedirs_shared} ${sourcedirs_2D})
endif(BUILD_WITH_IMAGE_PROCESS)


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

ADD_LIBRARY(sphinxsys_2d SHARED ${SCR_FILES})
ADD_LIBRARY(sphinxsys_static_2d STATIC ${SCR_FILES})

SET_TARGET_PROPERTIES(sphinxsys_static_2d PROPERTIES OUTPUT_NAME "sphinxsys_2d")
#SET_TARGET_PROPERTIES(sphinxsys PROPERTIES VERSION 1.0 SOVERSION 0)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

if(MSVC)
    target_link_libraries(sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES})
else(MSVC)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    	target_link_libraries(sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} ${Boost_LIBRARIES} stdc++)
	else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		target_link_libraries(sphinxsys_2d ${TBB_LIBRARYS} ${Simbody_LIBRARIES} stdc++ stdc++fs)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(MSVC)

INSTALL(TARGETS sphinxsys_2d sphinxsys_static_2d
RUNTIME DESTINATION 2d_code/bin
LIBRARY DESTINATION 2d_code/lib
ARCHIVE DESTINATION 2d_code/lib)

FILE(GLOB_RECURSE hpp_headers  ${PROJECT_SOURCE_DIR}/src/shared/*.hpp ${PROJECT_SOURCE_DIR}/src/for_2D_build/*.hpp)
INSTALL(FILES ${hpp_headers} DESTINATION 2d_code/include)
