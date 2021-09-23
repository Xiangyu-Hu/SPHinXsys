FILE(GLOB_RECURSE SPHINXsysHeaderPathListShared ${SPHINXSYS_PROJECT_DIR}/SPHINXsys/src/shared/*.h)
FILE(GLOB_RECURSE SPHINXsysHeaderPathListFor2DBuild ${SPHINXSYS_PROJECT_DIR}/SPHINXsys/src/for_2D_build/*.h)
FILE(GLOB_RECURSE SPHINXsysHeaderPathListFor2DBuildHpp ${CMAKE_CURRENT_SOURCE_DIR}/SPHINXsys/src/for_2D_build/*.hpp)

SET(SPHINXsysHeaderPath "")
FOREACH(file_path ${SPHINXsysHeaderPathListShared} ${SPHINXsysHeaderPathListFor2DBuild} ${SPHINXsysHeaderPathListFor2DBuildHpp})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(SPHINXsysHeaderPath ${SPHINXsysHeaderPath} ${dir_path})
ENDFOREACH()

if(BUILD_WITH_SIMBODY)
    include("${SPHINXSYS_PROJECT_DIR}/SPHINXsys/cmake/Simbody_header_directories.cmake")
    FOREACH(simbody_header_path ${SIMBODY_HEADER_DIRECTORIES})
        SET(SPHINXsysHeaderPath ${SPHINXsysHeaderPath} ${simbody_header_path})
    ENDFOREACH()
endif()

LIST(REMOVE_DUPLICATES SPHINXsysHeaderPath)
#message(STATUS ${SPHINXsysHeaderPath})

INCLUDE_DIRECTORIES("${SPHINXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/bin")