if(BUILD_WITH_SIMBODY)
    include(Simbody_header_directories)
endif()

MACRO(HEADER_DIRECTORIES_SHARED return_list)
    FILE(GLOB_RECURSE new_list  ${PROJECT_SOURCE_DIR}/src/shared/*.h)
    FILE(GLOB_RECURSE new_list_hpp  ${PROJECT_SOURCE_DIR}/src/shared/*.hpp)
   SET(dir_list "")
    FOREACH(file_path ${new_list} ${new_list_hpp})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)

    if(BUILD_WITH_SIMBODY)
        FOREACH(simbody_header_path ${SIMBODY_HEADER_DIRECTORIES})
            SET(dir_list ${dir_list} ${simbody_header_path})
        ENDFOREACH()
    endif()
    
    SET(${return_list} ${dir_list})
ENDMACRO()

MACRO(SOURCE_DIRECTORIES_SHARED return_list)
    FILE(GLOB_RECURSE  new_list  ${PROJECT_SOURCE_DIR}/src/shared/*.cpp)
    FILE(GLOB_RECURSE  new_list_c  ${PROJECT_SOURCE_DIR}/src/shared/*.c)
    SET(dir_list "")
    FOREACH(file_path ${new_list} ${new_list_c})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)

        if(BUILD_WITH_SIMBODY)
            SET(dir_list ${dir_list} ${dir_path})
        else() # if Simbody is not built, ignore the simbody and clapack folders
            if(NOT dir_list MATCHES ("/simbody/+" OR "/clapack_for_SPHinXsys/+") )
                SET(dir_list ${dir_list} ${dir_path})
            endif()
        endif()

    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

