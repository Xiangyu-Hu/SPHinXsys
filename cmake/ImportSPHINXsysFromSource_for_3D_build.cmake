FILE(GLOB_RECURSE SPHINXsysHeaderPathListShared ${CMAKE_SOURCE_DIR}/SPHINXsys/src/shared/*.h)
FILE(GLOB_RECURSE SPHINXsysHeaderPathListFor3DBuild ${CMAKE_SOURCE_DIR}/SPHINXsys/src/for_3D_build/*.h)
FILE(GLOB_RECURSE SPHINXsysHeaderPathListFor3DBuildHpp ${CMAKE_SOURCE_DIR}/SPHINXsys/src/for_3D_build/*.hpp)

SET(SPHINXsysHeaderPath "")
FOREACH(file_path ${SPHINXsysHeaderPathListShared} ${SPHINXsysHeaderPathListFor3DBuild} ${SPHINXsysHeaderPathListFor3DBuildHpp})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(SPHINXsysHeaderPath ${SPHINXsysHeaderPath} ${dir_path})
ENDFOREACH()
LIST(REMOVE_DUPLICATES SPHINXsysHeaderPath)

INCLUDE_DIRECTORIES("${SPHINXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHINXsys/bin")