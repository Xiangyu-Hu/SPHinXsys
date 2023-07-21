add_library(opencascade_modules aortic_valve_simulation/relax_dynamics_surface.cpp aortic_valve_simulation/relax_dynamics_surface.h aortic_valve_simulation/surface_shape.cpp aortic_valve_simulation/surface_shape.h aortic_valve_simulation/vector.cpp aortic_valve_simulation/vector.h)
target_include_directories(opencascade_modules PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/aortic_valve_simulation> $<INSTALL_INTERFACE:include/aortic_valve_simulation>)
target_link_libraries(opencascade_modules PUBLIC sphinxsys_3d)

# ## OCCT
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
find_package(OpenCASCADE CONFIG REQUIRED)
target_link_libraries(sphinxsys_core INTERFACE TKSTEP TKSTEPAttr TKSTEP209 TKSTEPBase TKXSBase TKShHealing TKTopAlgo TKGeomAlgo TKBRep TKGeomBase TKG3d TKG2d TKMath TKernel ) 
endif()