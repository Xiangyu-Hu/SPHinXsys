FILE(GLOB_RECURSE SPHINXsysHeaderPathList ${SPHINXSYS_PROJECT_DIR}/SPHINXsys/src/*.h)

SET(SPHINXsysHeaderPath "")
FOREACH(file_path ${SPHINXsysHeaderPathList})
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

INCLUDE_DIRECTORIES("${SPHINXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/bin")