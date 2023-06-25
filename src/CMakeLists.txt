file(GLOB_RECURSE SPHINXSYS_SHARED_HEADERS CONFIGURE_DEPENDS shared/*.h shared/*.hpp)
file(GLOB_RECURSE SPHINXSYS_SHARED_SOURCES CONFIGURE_DEPENDS shared/*.cpp)

if(SPHINXSYS_2D)
    add_library(sphinxsys_2d)
    #Make sure in-tree projects can reference this as SPHinXsys::sphinxsys_2d
    #to match the installed target names
    add_library(SPHinXsys::sphinxsys_2d ALIAS sphinxsys_2d)

    file(GLOB_RECURSE SPHINXSYS_2D_HEADERS CONFIGURE_DEPENDS for_2D_build/*.h for_2D_build/*.hpp)
    target_sources(sphinxsys_2d PUBLIC ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_2D_HEADERS})

    foreach(FILEPATH ${SPHINXSYS_SHARED_HEADERS})
        get_filename_component(DIR ${FILEPATH} DIRECTORY)
        list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
    endforeach()
    foreach(FILEPATH ${SPHINXSYS_2D_HEADERS})
        get_filename_component(DIR ${FILEPATH} DIRECTORY)
        list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
    endforeach()
    list(REMOVE_DUPLICATES SPHINXSYS_2D_INCLUDE_DIRS)
    foreach(DIR ${SPHINXSYS_2D_INCLUDE_DIRS})
        target_include_directories(sphinxsys_2d PUBLIC $<BUILD_INTERFACE:${DIR}>)
    endforeach()
    target_include_directories(sphinxsys_2d PUBLIC $<INSTALL_INTERFACE:include/sphinxsys_2d>)

    file(GLOB_RECURSE SPHINXSYS_2D_SOURCES CONFIGURE_DEPENDS for_2D_build/*.cpp)
    target_sources(sphinxsys_2d PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_2D_SOURCES})

    target_link_libraries(sphinxsys_2d PUBLIC sphinxsys_core)
endif()

if(SPHINXSYS_3D)

    add_library(sphinxsys_3d)
    #Make sure in-tree projects can reference this as SPHinXsys::sphinxsys_3d
    #to match the installed target names
    add_library(SPHinXsys::sphinxsys_3d ALIAS sphinxsys_3d)

    file(GLOB_RECURSE SPHINXSYS_3D_HEADERS CONFIGURE_DEPENDS for_3D_build/*.h for_3D_build/*.hpp)
    target_sources(sphinxsys_3d PUBLIC ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_3D_HEADERS})


    foreach(FILEPATH ${SPHINXSYS_SHARED_HEADERS})
      get_filename_component(DIR ${FILEPATH} DIRECTORY)
      list(APPEND SPHINXSYS_3D_INCLUDE_DIRS ${DIR})
    endforeach()
    foreach(FILEPATH ${SPHINXSYS_3D_HEADERS})
      get_filename_component(DIR ${FILEPATH} DIRECTORY)
      list(APPEND SPHINXSYS_3D_INCLUDE_DIRS ${DIR})
    endforeach()
    list(REMOVE_DUPLICATES SPHINXSYS_3D_INCLUDE_DIRS)
    foreach(DIR ${SPHINXSYS_3D_INCLUDE_DIRS})
        target_include_directories(sphinxsys_3d PUBLIC $<BUILD_INTERFACE:${DIR}>)
    endforeach()
    target_include_directories(sphinxsys_3d PUBLIC $<INSTALL_INTERFACE:include/sphinxsys_3d>)

    file(GLOB_RECURSE SPHINXSYS_3D_SOURCES CONFIGURE_DEPENDS for_3D_build/*.cpp)
    target_sources(sphinxsys_3d PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_3D_SOURCES})

    target_link_libraries(sphinxsys_3d PUBLIC sphinxsys_core)
endif()
