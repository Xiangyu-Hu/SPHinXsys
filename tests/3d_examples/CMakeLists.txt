SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT EMSCRIPTEN)
    FOREACH(subdir ${SUBDIRS})
            ADD_SUBDIRECTORY(${subdir})
    ENDFOREACH()
else()
    add_subdirectory(test_3d_passive_cantilever_neohookean)
endif()