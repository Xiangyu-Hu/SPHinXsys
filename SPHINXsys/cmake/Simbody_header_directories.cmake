set(SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS "")
if(NOT WASM_BUILD)
    set(SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS "${CMAKE_SOURCE_DIR}") #SPHinXsys
else()
    set(SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS "${CMAKE_SOURCE_DIR}/SPHinXsys_Virtonomy/SPHinXsys") #SPH_to_JS/SPHinXsys_Virtonomy/SPHinXsys
endif()

SET(SIMBODY_HEADER_DIRECTORIES
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/BigMatrix/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/Geometry/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/Mechanics/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/Polynomial/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/Random/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/Scalar/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/Simulation/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKcommon/SmallMatrix/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKmath/Geometry/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/Simbody_3_7/SimTKmath/Integrators/include"
)