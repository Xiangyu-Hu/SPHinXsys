FILE(GLOB_RECURSE SPHinXsysHeaderPathList ${SPHINXSYS_PROJECT_DIR}/src/*.h)

SET(SPHinXsysHeaderPath "")

FOREACH(file_path ${SPHinXsysHeaderPathList})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(SPHinXsysHeaderPath ${SPHinXsysHeaderPath} ${dir_path})
ENDFOREACH()

LIST(REMOVE_DUPLICATES SPHinXsysHeaderPath)

INCLUDE_DIRECTORIES("${SPHinXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHinXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHinXsys/bin")