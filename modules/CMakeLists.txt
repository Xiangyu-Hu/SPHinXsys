include(ImportSPHINXsysFromSource_for_3D_build)

if(NOT STATIC_BUILD)
    add_library(structural_simulation_module SHARED structural_simulation/structural_simulation_class.cpp)
    target_link_libraries(structural_simulation_module sphinxsys_3d)
else()
    add_library(structural_simulation_module STATIC structural_simulation/structural_simulation_class.cpp)
    target_link_libraries(structural_simulation_module sphinxsys_static_3d)
endif()

target_include_directories(structural_simulation_module PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/structural_simulation ${SPHINXsysHeaderPath})

if(BUILD_WITH_SIMBODY)
    include("${SPHINXSYS_PROJECT_DIR}/SPHINXsys/cmake/Simbody_header_directories.cmake")
    target_include_directories(structural_simulation_module PUBLIC ${SIMBODY_HEADER_DIRECTORIES})
endif()


