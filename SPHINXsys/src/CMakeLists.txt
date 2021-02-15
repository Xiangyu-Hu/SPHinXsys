# shared build 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake) # main (top) cmake dir
include(Dirsearch_shared)

## prepare dirctoriesfor head and source files
HEADER_DIRECTORIES_SHARED(headdirs_shared)
SOURCE_DIRECTORIES_SHARED(sourcedirs_shared)

SET(usefuldirs ${headdirs_shared} ${sourcedirs_shared})
LIST(REMOVE_DUPLICATES usefuldirs)

SET(usefulsubdirs ${usefuldirs})
LIST(REMOVE_ITEM usefulsubdirs ${CMAKE_CURRENT_SOURCE_DIR})

#Add all useful subdirectories
FOREACH(subdir_path ${usefulsubdirs})
	#message(STATUS ${subdir_path})
	ADD_SUBDIRECTORY(${subdir_path})
ENDFOREACH()

#ADD_SUBDIRECTORY(for_2D_build)
ADD_SUBDIRECTORY(for_3D_build)
ADD_SUBDIRECTORY(shared)
ADD_SUBDIRECTORY(oneTBB)

include_directories(oneTBB)