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

if(NOT ONLY_3D)
	ADD_SUBDIRECTORY(for_2D_build)
endif()
ADD_SUBDIRECTORY(for_3D_build)
ADD_SUBDIRECTORY(shared)

if(BUILD_WITH_IMAGE_PROCESS)
add_subdirectory(image_processing)
endif()