set(SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS "")
if(NOT WASM_BUILD)
    set(SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS "${CMAKE_SOURCE_DIR}") #SPHinXsys
else()
    set(SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS "${CMAKE_SOURCE_DIR}/SPHinXsys_Virtonomy/SPHinXsys") #SPH_to_JS/SPHinXsys_Virtonomy/SPHinXsys
endif()

SET(SIMBODY_HEADER_DIRECTORIES
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/BigMatrix/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/Geometry/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/Mechanics/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/Polynomial/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/Random/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/Scalar/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/Simulation/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKcommon/SmallMatrix/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKmath/Geometry/include"
    "${SHP_SOURCE_DIRECTORY_FOR_INCLUDE_PATHS}/SPHINXsys/src/shared/simbody/SimTKmath/Integrators/include"
)