if(NOT SPHINXSYS_MODULE_OPENCASCADE OR ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    return()
endif()

find_package(OpenCASCADE CONFIG REQUIRED)
# Must be shared library if OpenCASCADE is built static otherwise coplyeft OpenCASCADE license propagates
get_target_property(OCCT_LIBRARY_TYPE TKernel TYPE) # https://cmake.org/cmake/help/latest/prop_tgt/TYPE.html
if(OCCT_LIBRARY_TYPE MATCHES "STATIC_LIBRARY")
    set(SPHINXSYS_OPENCASCADE_LIBRARY_TYPE SHARED)
endif()
add_library(sphinxsys_opencascade
    ${SPHINXSYS_OPENCASCADE_LIBRARY_TYPE}
    opencascade/relax_dynamics_surface.cpp 
    opencascade/relax_dynamics_surface.h 
    opencascade/surface_shape.cpp 
    opencascade/surface_shape.h 
    opencascade/vector.cpp 
    opencascade/vector.h
)
target_include_directories(sphinxsys_opencascade PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/opencascade>)
target_link_libraries(sphinxsys_opencascade PUBLIC sphinxsys_3d)

# Must list all dependent packages given here https://dev.opencascade.org/doc/refman/html/toolkit_tkstep.html because improper CMake target link visibility
target_link_libraries(sphinxsys_opencascade PRIVATE TKDESTEP TKXSBase TKShHealing TKTopAlgo TKGeomAlgo TKBRep TKGeomBase) 
target_link_libraries(sphinxsys_opencascade PUBLIC TKG3d TKG2d TKMath TKernel) 

if(SPHINXSYS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()