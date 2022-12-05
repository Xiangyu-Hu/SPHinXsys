set(simbody_root "${SPHINXSYS_PROJECT_DIR}/3rd_party/simbody")

set(SIMBODY_HEADER_DIRECTORIES
    "${simbody_root}/SimTKcommon/BigMatrix/include"
    "${simbody_root}/SimTKcommon/Geometry/include"
    "${simbody_root}/SimTKcommon/Mechanics/include"
    "${simbody_root}/SimTKcommon/Polynomial/include"
    "${simbody_root}/SimTKcommon/Random/include"
    "${simbody_root}/SimTKcommon/Scalar/include"
    "${simbody_root}/SimTKcommon/Simulation/include"
    "${simbody_root}/SimTKcommon/SmallMatrix/include"
    "${simbody_root}/SimTKcommon/include"

    "${simbody_root}/SimTKmath/Geometry/include"
    "${simbody_root}/SimTKmath/Integrators/include"
    "${simbody_root}/SimTKmath/include"

    "${simbody_root}/Simbody/include"
    "${simbody_root}/Simbody/include/simbody/internal"
    "${simbody_root}/Simbody/Visualizer/include"
)