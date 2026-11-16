file(GLOB_RECURSE SPHINXSYS_SHARED_HEADERS CONFIGURE_DEPENDS shared/*.h shared/*.hpp  shared/*.hxx)
file(GLOB_RECURSE SPHINXSYS_SHARED_SOURCES CONFIGURE_DEPENDS shared/*.cpp)

if(SPHINXSYS_USE_SYCL)
file(GLOB_RECURSE SPHINXSYS_SYCL_HEADERS CONFIGURE_DEPENDS src_sycl/*.h src_sycl/*.hpp)
file(GLOB_RECURSE SPHINXSYS_SYCL_SOURCES CONFIGURE_DEPENDS src_sycl/*.cpp)
endif()

file(GLOB_RECURSE SIMBODY_SYSTEM_SOURCES CONFIGURE_DEPENDS simbody_system/*.cpp) 

if(SPHINXSYS_2D)
    add_library(sphinxsys_2d)
    #Make sure in-tree projects can reference this as SPHinXsys::sphinxsys_2d
    #to match the installed target names
    add_library(SPHinXsys::sphinxsys_2d ALIAS sphinxsys_2d)

    file(GLOB_RECURSE SPHINXSYS_2D_HEADERS CONFIGURE_DEPENDS for_2D_build/*.h for_2D_build/*.hpp)
    file(GLOB_RECURSE SPHINXSYS_2D_SOURCES CONFIGURE_DEPENDS for_2D_build/*.cpp)
    
    foreach(FILEPATH ${SPHINXSYS_SHARED_HEADERS})
        get_filename_component(DIR ${FILEPATH} DIRECTORY)
        list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
    endforeach()
    
    foreach(FILEPATH ${SPHINXSYS_2D_HEADERS})
        get_filename_component(DIR ${FILEPATH} DIRECTORY)
        list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
    endforeach()  
    
    if(SPHINXSYS_USE_SYCL)
        foreach(FILEPATH ${SPHINXSYS_SYCL_HEADERS})
            get_filename_component(DIR ${FILEPATH} DIRECTORY)
            list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
        endforeach()

        target_sources(sphinxsys_2d PUBLIC ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_2D_HEADERS} ${SPHINXSYS_SYCL_HEADERS})
        target_sources(sphinxsys_2d PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_2D_SOURCES} ${SPHINXSYS_SYCL_SOURCES})   
    else()
        target_sources(sphinxsys_2d PUBLIC ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_2D_HEADERS})
        target_sources(sphinxsys_2d PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_2D_SOURCES})
    endif()


    list(REMOVE_DUPLICATES SPHINXSYS_2D_INCLUDE_DIRS)
 
    add_library(simbody_system_2d OBJECT ${SIMBODY_SYSTEM_SOURCES})

    foreach(DIR ${SPHINXSYS_2D_INCLUDE_DIRS})
        target_include_directories(simbody_system_2d PRIVATE ${DIR})
    endforeach()

    if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        target_compile_options(simbody_system_2d PRIVATE -fp-model=precise -fno-fast-math)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(simbody_system_2d PRIVATE)
    endif()
    find_package(Eigen3 CONFIG REQUIRED)
    target_link_libraries(simbody_system_2d PRIVATE sphinxsys_core)

    foreach(DIR ${SPHINXSYS_2D_INCLUDE_DIRS})
        target_include_directories(sphinxsys_2d PUBLIC $<BUILD_INTERFACE:${DIR}>)
    endforeach()

    target_include_directories(sphinxsys_2d PUBLIC $<INSTALL_INTERFACE:include/sphinxsys_2d>)
    
    target_link_libraries(sphinxsys_2d PUBLIC sphinxsys_lib)
    target_sources(sphinxsys_2d PRIVATE $<TARGET_OBJECTS:simbody_system_2d>)
    target_compile_definitions(sphinxsys_2d PUBLIC SPH_NDIM=2)
endif()

if(SPHINXSYS_3D)

    add_library(sphinxsys_3d)
    #Make sure in-tree projects can reference this as SPHinXsys::sphinxsys_3d
    #to match the installed target names
    add_library(SPHinXsys::sphinxsys_3d ALIAS sphinxsys_3d)

    file(GLOB_RECURSE SPHINXSYS_3D_HEADERS CONFIGURE_DEPENDS for_3D_build/*.h for_3D_build/*.hpp)
    file(GLOB_RECURSE SPHINXSYS_3D_SOURCES CONFIGURE_DEPENDS for_3D_build/*.cpp)

    foreach(FILEPATH ${SPHINXSYS_SHARED_HEADERS})
        get_filename_component(DIR ${FILEPATH} DIRECTORY)
        list(APPEND SPHINXSYS_3D_INCLUDE_DIRS ${DIR})
    endforeach()

    foreach(FILEPATH ${SPHINXSYS_3D_HEADERS})
        get_filename_component(DIR ${FILEPATH} DIRECTORY)
        list(APPEND SPHINXSYS_3D_INCLUDE_DIRS ${DIR})
    endforeach()

    if(SPHINXSYS_USE_SYCL)
        foreach(FILEPATH ${SPHINXSYS_SYCL_HEADERS})
            get_filename_component(DIR ${FILEPATH} DIRECTORY)
            list(APPEND SPHINXSYS_3D_INCLUDE_DIRS ${DIR})
        endforeach()

        target_sources(sphinxsys_3d PUBLIC ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_3D_HEADERS} ${SPHINXSYS_SYCL_HEADERS})
        target_sources(sphinxsys_3d PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_3D_SOURCES} ${SPHINXSYS_SYCL_SOURCES})
    else()
    target_sources(sphinxsys_3d PUBLIC ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_3D_HEADERS})
    target_sources(sphinxsys_3d PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_3D_SOURCES})
    endif()

    list(REMOVE_DUPLICATES SPHINXSYS_3D_INCLUDE_DIRS)

    add_library(simbody_system_3d OBJECT ${SIMBODY_SYSTEM_SOURCES})

    foreach(DIR ${SPHINXSYS_3D_INCLUDE_DIRS})
        target_include_directories(simbody_system_3d PRIVATE ${DIR})
    endforeach()


    if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        target_compile_options(simbody_system_3d PRIVATE -fp-model=precise -fno-fast-math)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(simbody_system_3d PRIVATE)
    endif()
    target_link_libraries(simbody_system_3d PRIVATE sphinxsys_core)


    foreach(DIR ${SPHINXSYS_3D_INCLUDE_DIRS})
        target_include_directories(sphinxsys_3d PUBLIC $<BUILD_INTERFACE:${DIR}>)
    endforeach()
    target_include_directories(sphinxsys_3d PUBLIC $<INSTALL_INTERFACE:include/sphinxsys_3d>)

    target_link_libraries(sphinxsys_3d PUBLIC sphinxsys_lib)
    target_sources(sphinxsys_3d PRIVATE $<TARGET_OBJECTS:simbody_system_3d>)
    target_compile_definitions(sphinxsys_3d PUBLIC SPH_NDIM=3)
endif()
