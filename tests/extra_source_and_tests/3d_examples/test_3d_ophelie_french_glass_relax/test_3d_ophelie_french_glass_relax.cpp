/**
 * @file test_3d_ophelie_french_glass_relax.cpp
 * @brief French reduced glass: TriangleMeshShapeCylinder + SphinxSys SYCL-CK particle relaxation.
 *
 * Writes Reload.xml (GlassBody) for test_3d_ophelie_french_reduced EM runs.
 * Geometry matches Jacoutot-inspired reduced case (D=650 mm cylinder, mesh resolution configurable).
 */
#include "electromagnetic_ophelie_french_glass_mesh_relax.h"
#include "electromagnetic_ophelie_progress.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;

int main(int ac, char *av[])
{
    logOphelieRunContext();
    std::cout << "[ophelie] French glass mesh relax (TriangleMeshShapeCylinder + SYCL-CK)\n"
              << "[ophelie] EM step: test_3d_ophelie_french_reduced --reload=1 ...\n";

    OphelieFrenchReducedCaseParams french;
    OphelieParameters params;
    applyFrenchReducedDefaults(params, french);

    const StdVec<std::string> filtered = filterFrenchReducedCommandLine(ac, av, french);
    refreshFrenchReducedCoilStack(french);

    size_t relaxation_steps = 1000;
    size_t vtp_every = 100;
    StdVec<std::string> sph_arguments;
    sph_arguments.emplace_back(filtered.front());
    for (size_t i = 1; i < filtered.size(); ++i)
    {
        const std::string &arg = filtered[i];
        if (arg.rfind("--relax-steps=", 0) == 0)
        {
            relaxation_steps = static_cast<size_t>(std::atoi(arg.c_str() + 14));
            continue;
        }
        if (arg.rfind("--relax-vtp-every=", 0) == 0)
        {
            vtp_every = static_cast<size_t>(std::atoi(arg.c_str() + 18));
            continue;
        }
        sph_arguments.emplace_back(arg);
    }

    StdVec<char *> sph_av;
    sph_av.reserve(sph_arguments.size());
    for (auto &argument : sph_arguments)
    {
        sph_av.push_back(const_cast<char *>(argument.c_str()));
    }

    const Real pad = 3.0 * french.dp;
    const BoundingBoxd system_bounds = frenchGlassMeshRelaxDomainBounds(french, pad);

    IO::getEnvironment().resetReloadFolder("./reload", true);

    SPHSystem sph_system(system_bounds, french.dp);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());

    printFrenchReducedCaseSummary(french);
    logFrenchReducedCoilGeometry(french);

    SolidBody glass_body(
        sph_system, makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                       french.glass_radius, french.glass_half_height,
                                                                       french.glass_mesh_resolution));
    glass_body.defineAdaptation<SPHAdaptation>(1.0, 1.0);
    glass_body.defineMatterMaterial<Solid>();

#if SPHINXSYS_USE_SYCL
    LevelSetShape &level_set_shape = defineFrenchGlassMeshLevelSet(glass_body, true);
#else
    std::cout << "[ophelie] ERROR: SYCL required for mesh-cylinder relax." << std::endl;
    return 1;
#endif

    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    const size_t n0 = glass_body.getBaseParticles().TotalRealParticles();
    std::cout << "[ophelie] lattice particles n=" << n0 << " dp=" << french.dp << std::endl;

    runFrenchGlassSphinxsysStyleRelax(sph_system, glass_body, level_set_shape, relaxation_steps, vtp_every);

    ReloadParticleIO write_reload(glass_body);
    write_reload.writeToFile(0);

    const std::string output_folder = IO::getEnvironment().OutputFolder();
    logOphelieOutputArtifact(output_folder + "/GlassBody_ite_0000001000.vtp");
    logOphelieOutputArtifact(IO::getEnvironment().ReloadFolder() + "/Reload.xml");

    std::cout << "test_3d_ophelie_french_glass_relax finished. n=" << glass_body.getBaseParticles().TotalRealParticles()
              << " Reload.xml (GlassBody) -> " << IO::getEnvironment().ReloadFolder()
              << " | relax VTP snapshots -> " << output_folder << "/GlassBody_ite_*.vtp" << std::endl;
    return 0;
}
