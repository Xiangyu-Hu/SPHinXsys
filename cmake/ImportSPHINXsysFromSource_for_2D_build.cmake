FILE(GLOB_RECURSE SPHINXsysHeaderPathListShared ${CMAKE_SOURCE_DIR}/SPHINXsys/src/shared/*.h)
FILE(GLOB_RECURSE SPHINXsysHeaderPathListFor2DBuild ${CMAKE_SOURCE_DIR}/SPHINXsys/src/for_2D_build/*.h)
FILE(GLOB_RECURSE SPHINXsysHeaderPathListFor2DBuildHpp ${CMAKE_SOURCE_DIR}/SPHINXsys/src/for_2D_build/*.hpp)

include("${CMAKE_SOURCE_DIR}/SPHINXsys/cmake/Simbody_header_directories.cmake")

SET(SPHINXsysHeaderPath "")
FOREACH(file_path ${SPHINXsysHeaderPathListShared} ${SPHINXsysHeaderPathListFor2DBuild} ${SPHINXsysHeaderPathListFor2DBuildHpp})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(SPHINXsysHeaderPath ${SPHINXsysHeaderPath} ${dir_path})
ENDFOREACH()

FOREACH(simbody_header_path ${SIMBODY_HEADER_DIRECTORIES})
    SET(SPHINXsysHeaderPath ${SPHINXsysHeaderPath} ${simbody_header_path})
ENDFOREACH()

if(NOT oneTBB_BUILD MATCHES "off")
    include("${CMAKE_SOURCE_DIR}/SPHINXsys/cmake/oneTBB_header_directory.cmake")
    SET(SPHINXsysHeaderPath ${SPHINXsysHeaderPath} ${ONETBB_HEADER_DIRECTORY})
endif(NOT oneTBB_BUILD MATCHES "off")

LIST(REMOVE_DUPLICATES SPHINXsysHeaderPath)
#message(STATUS ${SPHINXsysHeaderPath})

INCLUDE_DIRECTORIES("${SPHINXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/bin")