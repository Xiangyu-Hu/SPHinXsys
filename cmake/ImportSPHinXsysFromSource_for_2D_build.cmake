FILE(GLOB_RECURSE SPHinXsysHeaderPathListShared ${SPHINXSYS_PROJECT_DIR}/src/shared/*.h)
FILE(GLOB_RECURSE SPHinXsysHeaderPathListFor2DBuild ${SPHINXSYS_PROJECT_DIR}/src/for_2D_build/*.h)
FILE(GLOB_RECURSE SPHinXsysHeaderPathListFor2DBuildHpp ${CMAKE_CURRENT_SOURCE_DIR}/src/for_2D_build/*.hpp)

SET(SPHinXsysHeaderPath "")

FOREACH(file_path ${SPHinXsysHeaderPathListShared} ${SPHinXsysHeaderPathListFor2DBuild} ${SPHinXsysHeaderPathListFor2DBuildHpp})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(SPHinXsysHeaderPath ${SPHinXsysHeaderPath} ${dir_path})
ENDFOREACH()

LIST(REMOVE_DUPLICATES SPHinXsysHeaderPath)

# message(STATUS ${SPHinXsysHeaderPath})
INCLUDE_DIRECTORIES("${SPHinXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHinXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHinXsys/bin")