# GPU-paper branch: only ship recording helpers under shared/gpu_paper.
# Legacy extras (turbulence UDFs etc.) break under icpx (esp. FP32) and are not SYCL-safe.
# Set SPHINXSYS_BUILD_LEGACY_EXTRA_SOURCES=ON to restore the full extra_src tree.

if(SPHINXSYS_USE_SYCL OR NOT SPHINXSYS_BUILD_LEGACY_EXTRA_SOURCES)
    file(GLOB_RECURSE GPU_PAPER_HEADERS CONFIGURE_DEPENDS shared/gpu_paper/*.h shared/gpu_paper/*.hpp)
    file(GLOB_RECURSE GPU_PAPER_SOURCES CONFIGURE_DEPENDS shared/gpu_paper/*.cpp)

    if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)
        add_library(extra_sources_2d)
        target_sources(extra_sources_2d PUBLIC ${GPU_PAPER_HEADERS})
        target_sources(extra_sources_2d PRIVATE ${GPU_PAPER_SOURCES})
        target_include_directories(extra_sources_2d PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/shared/gpu_paper>
            $<INSTALL_INTERFACE:include/extra_sources_2d>)
        target_link_libraries(extra_sources_2d PUBLIC sphinxsys_2d)
    endif()

    if(SPHINXSYS_3D AND SPHINXSYS_BUILD_3D_EXAMPLES)
        add_library(extra_sources_3d)
        target_sources(extra_sources_3d PUBLIC ${GPU_PAPER_HEADERS})
        target_sources(extra_sources_3d PRIVATE ${GPU_PAPER_SOURCES})
        target_include_directories(extra_sources_3d PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/shared/gpu_paper>
            $<INSTALL_INTERFACE:include/extra_sources_3d>)
        target_link_libraries(extra_sources_3d PUBLIC sphinxsys_3d)
    endif()
else()
    file(GLOB_RECURSE EXTRA_SOURCES_SHARED_HEADERS CONFIGURE_DEPENDS shared/*.h shared/*.hpp)
    file(GLOB_RECURSE EXTRA_SOURCES_SHARED_SOURCES CONFIGURE_DEPENDS shared/*.cpp)

    if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)
        add_library(extra_sources_2d)

        file(GLOB_RECURSE EXTRA_SOURCES_2D_HEADERS CONFIGURE_DEPENDS for_2D_build/*.h for_2D_build/*.hpp)
        target_sources(extra_sources_2d PUBLIC ${EXTRA_SOURCES_SHARED_HEADERS} ${EXTRA_SOURCES_2D_HEADERS})

        foreach(FILEPATH ${EXTRA_SOURCES_SHARED_HEADERS})
            get_filename_component(DIR ${FILEPATH} DIRECTORY)
            list(APPEND EXTRA_SOURCES_2D_INCLUDE_DIRS ${DIR})
        endforeach()
        foreach(FILEPATH ${EXTRA_SOURCES_2D_HEADERS})
            get_filename_component(DIR ${FILEPATH} DIRECTORY)
            list(APPEND EXTRA_SOURCES_2D_INCLUDE_DIRS ${DIR})
        endforeach()
        list(REMOVE_DUPLICATES EXTRA_SOURCES_2D_INCLUDE_DIRS)
        foreach(DIR ${EXTRA_SOURCES_2D_INCLUDE_DIRS})
            target_include_directories(extra_sources_2d PUBLIC $<BUILD_INTERFACE:${DIR}>)
        endforeach()
        target_include_directories(extra_sources_2d PUBLIC $<INSTALL_INTERFACE:include/extra_sources_2d>)

        file(GLOB_RECURSE EXTRA_SOURCES_2D_SOURCES CONFIGURE_DEPENDS for_2D_build/*.cpp)
        target_sources(extra_sources_2d PRIVATE ${EXTRA_SOURCES_SHARED_SOURCES} ${EXTRA_SOURCES_2D_SOURCES})

        target_link_libraries(extra_sources_2d PUBLIC sphinxsys_2d)
    endif()

    if(SPHINXSYS_3D AND SPHINXSYS_BUILD_3D_EXAMPLES)
        add_library(extra_sources_3d)

        file(GLOB_RECURSE EXTRA_SOURCES_3D_HEADERS CONFIGURE_DEPENDS for_3D_build/*.h for_3D_build/*.hpp)
        target_sources(extra_sources_3d PUBLIC ${EXTRA_SOURCES_SHARED_HEADERS} ${EXTRA_SOURCES_3D_HEADERS})

        foreach(FILEPATH ${EXTRA_SOURCES_SHARED_HEADERS})
            get_filename_component(DIR ${FILEPATH} DIRECTORY)
            list(APPEND EXTRA_SOURCES_3D_INCLUDE_DIRS ${DIR})
        endforeach()
        foreach(FILEPATH ${EXTRA_SOURCES_3D_HEADERS})
            get_filename_component(DIR ${FILEPATH} DIRECTORY)
            list(APPEND EXTRA_SOURCES_3D_INCLUDE_DIRS ${DIR})
        endforeach()
        list(REMOVE_DUPLICATES EXTRA_SOURCES_3D_INCLUDE_DIRS)
        foreach(DIR ${EXTRA_SOURCES_3D_INCLUDE_DIRS})
            target_include_directories(extra_sources_3d PUBLIC $<BUILD_INTERFACE:${DIR}>)
        endforeach()
        target_include_directories(extra_sources_3d PUBLIC $<INSTALL_INTERFACE:include/extra_sources_3d>)

        file(GLOB_RECURSE EXTRA_SOURCES_3D_SOURCES CONFIGURE_DEPENDS for_3D_build/*.cpp)
        target_sources(extra_sources_3d PRIVATE ${EXTRA_SOURCES_SHARED_SOURCES} ${EXTRA_SOURCES_3D_SOURCES})

        target_link_libraries(extra_sources_3d PUBLIC sphinxsys_3d)
    endif()
endif()
