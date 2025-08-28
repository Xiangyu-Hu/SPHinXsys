FILE(GLOB_RECURSE SPHinXsysHeaderPathListShared ${SPHINXSYS_PROJECT_DIR}/src/shared/*.h)
FILE(GLOB_RECURSE SPHinXsysHeaderPathListFor3DBuild ${SPHINXSYS_PROJECT_DIR}/src/for_3D_build/*.h)
FILE(GLOB_RECURSE SPHinXsysHeaderPathListFor3DBuildHpp ${SPHINXSYS_PROJECT_DIR}/src/for_3D_build/*.hpp)

SET(SPHinXsysHeaderPath "")

FOREACH(file_path ${SPHinXsysHeaderPathListShared} ${SPHinXsysHeaderPathListFor3DBuild} ${SPHinXsysHeaderPathListFor3DBuildHpp})
    GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
    SET(SPHinXsysHeaderPath ${SPHinXsysHeaderPath} ${dir_path})
ENDFOREACH()

LIST(REMOVE_DUPLICATES SPHinXsysHeaderPath)

INCLUDE_DIRECTORIES("${SPHinXsysHeaderPath}")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHinXsys/lib")
LINK_DIRECTORIES("${CMAKE_BINARY_DIR}/SPHinXsys/bin")