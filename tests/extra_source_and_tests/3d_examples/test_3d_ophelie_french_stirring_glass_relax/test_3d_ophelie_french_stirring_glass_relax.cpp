/**
 * @file test_3d_ophelie_french_stirring_glass_relax.cpp
 * @brief French stirring CAD: SYCL level-set particle relaxation of paddle and/or molten glass.
 *
 * Structure follows sphinxsys-einsimo test_3d_stirring_relaxation_oil_two_phase:
 * per-body level set (refined for thin paddle blades) -> Lattice -> randomize ->
 * KernelGradientIntegral + LevelsetBounding relaxation -> reload.
 *
 * Bodies (select with --bodies=all or a comma-separated subset of glass,rotor,wall):
 *   Rotor        — stirring_paddle_2_z.stl
 *   GlassBody    — glass-z.stl minus the paddle
 *   WallBoundary — crucible shell around the melt, open at the free surface
 *
 * CAD units are mm, so --geometry-scale=0.001 (default) and no translation is applied.
 */
#include "electromagnetic_ophelie_french_stirring_geometry.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    (void)ac;
    (void)av;
    std::cerr << "test_3d_ophelie_french_stirring_glass_relax requires SYCL build (SPHINXSYS_USE_SYCL)." << std::endl;
    return 1;
#else
    OphelieFrenchStirringCaseParams stirring;
    applyFrenchReducedStirringDefaults(stirring);
    const StdVec<std::string> filtered = filterFrenchStirringCommandLine(ac, av, stirring);
    StdVec<char *> sph_av;
    for (const std::string &argument : filtered)
    {
        sph_av.push_back(const_cast<char *>(argument.c_str()));
    }

    const bool build_rotor = frenchStirringRelaxRotor(stirring);
    const bool build_glass = frenchStirringRelaxGlass(stirring);
    const bool build_wall = frenchStirringRelaxWall(stirring);
    if (!build_rotor && !build_glass && !build_wall)
    {
        std::cerr << "[ophelie][stirring-relax] --bodies must be \"all\" or a subset of glass,rotor,wall"
                  << std::endl;
        return 1;
    }

    const BoundingBoxd system_domain_bounds = frenchStirringRelaxDomainBounds(stirring);
    SPHSystem sph_system(system_domain_bounds, stirring.french.dp);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());
    IO::getEnvironment().resetReloadFolder("./reload", true);
    printFrenchStirringCaseSummary(stirring);
    std::cout << "[ophelie][stirring-relax] bodies=" << stirring.relax_bodies
              << " rotor_level_set_refinement=" << stirring.rotor_level_set_refinement
              << " glass_level_set_refinement=" << stirring.glass_level_set_refinement << std::endl;

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    StdVec<RealBody *> relax_bodies;
    StdVec<SPHBody *> reload_bodies;
    ParticleDynamicsGroup update_relation;
    ParticleDynamicsGroup relaxation_residual;
    ParticleDynamicsGroup surface_bounding;

    UniquePtr<SolidBody> rotor_keeper;
    UniquePtr<SolidBody> glass_keeper;
    UniquePtr<SolidBody> wall_keeper;
    UniquePtr<Inner<>> rotor_inner_keeper;
    UniquePtr<Inner<>> glass_inner_keeper;
    UniquePtr<Inner<>> wall_inner_keeper;
    UniquePtr<NearShapeSurface> rotor_near_surface_keeper;
    UniquePtr<NearShapeSurface> glass_near_surface_keeper;
    UniquePtr<NearShapeSurface> wall_near_surface_keeper;

    if (build_rotor)
    {
        rotor_keeper = makeUnique<SolidBody>(sph_system, makeShared<OphelieFrenchStirringRotorShape>("Rotor", stirring));
        SolidBody &rotor = *rotor_keeper;
        rotor.defineAdaptationRatios(1.3, 1.0);
        rotor.defineMatterMaterial<Solid>(stirring.rho_rotor);
        LevelSetShape &rotor_level_set_shape =
            rotor.defineBodyLevelSetShape(par_ck, stirring.rotor_level_set_refinement)
                .correctLevelSetSign()
                .cleanLevelSet();
        rotor.generateParticles<BaseParticles, Lattice>();

        rotor_inner_keeper = makeUnique<Inner<>>(rotor);
        rotor_near_surface_keeper = makeUnique<NearShapeSurface>(rotor);
        relax_bodies.push_back(&rotor);
        reload_bodies.push_back(&rotor);

        host_methods.addStateDynamics<RandomizeParticlePositionCK>(rotor).exec();
        update_relation.add(&main_methods.addRelationDynamics(*rotor_inner_keeper));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(*rotor_inner_keeper)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(rotor, rotor_level_set_shape));
        surface_bounding.add(&main_methods.addStateDynamics<LevelsetBounding>(*rotor_near_surface_keeper));
    }

    if (build_glass)
    {
        glass_keeper =
            makeUnique<SolidBody>(sph_system, makeShared<OphelieFrenchStirringGlassFluidShape>("GlassBody", stirring));
        SolidBody &glass = *glass_keeper;
        glass.defineAdaptation<SPHAdaptation>(1.0, 1.0); // same as test_3d_ophelie_french_glass_relax
        glass.defineMatterMaterial<Solid>();
        LevelSetShape &glass_level_set_shape =
            glass.defineBodyLevelSetShape(par_ck, stirring.glass_level_set_refinement)
                .correctLevelSetSign()
                .cleanLevelSet();
        glass.generateParticles<BaseParticles, Lattice>();

        glass_inner_keeper = makeUnique<Inner<>>(glass);
        glass_near_surface_keeper = makeUnique<NearShapeSurface>(glass);
        relax_bodies.push_back(&glass);
        reload_bodies.push_back(&glass);

        host_methods.addStateDynamics<RandomizeParticlePositionCK>(glass).exec();
        update_relation.add(&main_methods.addRelationDynamics(*glass_inner_keeper));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(*glass_inner_keeper)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(glass, glass_level_set_shape));
        surface_bounding.add(&main_methods.addStateDynamics<LevelsetBounding>(*glass_near_surface_keeper));
    }

    if (build_wall)
    {
        wall_keeper = makeUnique<SolidBody>(
            sph_system, makeShared<OphelieFrenchNaturalCrucibleWallShape>(
                            "WallBoundary", stirring.french, frenchStirringWallThickness(stirring),
                            stirring.wall_mesh_resolution));
        SolidBody &wall = *wall_keeper;
        wall.defineAdaptationRatios(1.3, 1.0);
        wall.defineMatterMaterial<Solid>();
        LevelSetShape &wall_level_set_shape =
            wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
        wall.generateParticles<BaseParticles, Lattice>();

        wall_inner_keeper = makeUnique<Inner<>>(wall);
        wall_near_surface_keeper = makeUnique<NearShapeSurface>(wall);
        relax_bodies.push_back(&wall);
        reload_bodies.push_back(&wall);

        host_methods.addStateDynamics<RandomizeParticlePositionCK>(wall).exec();
        update_relation.add(&main_methods.addRelationDynamics(*wall_inner_keeper));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(*wall_inner_keeper)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(wall, wall_level_set_shape));
        surface_bounding.add(&main_methods.addStateDynamics<LevelsetBounding>(*wall_near_surface_keeper));
    }

    ParticleDynamicsGroup update_configuration = main_methods.addCellLinkedListDynamics(relax_bodies) + update_relation;
    ReduceDynamicsGroup relaxation_scaling =
        main_methods.addReduceDynamics<ReduceMin, RelaxationScalingCK>(relax_bodies);
    ParticleDynamicsGroup update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(relax_bodies);
    update_particle_position = update_particle_position + surface_bounding;

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    auto &write_particle_reload_files = main_methods.addIODynamics<ReloadParticleIOCK>(reload_bodies);

    update_configuration.exec();
    body_state_recorder.writeToFile(0);

    std::cout << "[ophelie][stirring-relax] n_rotor="
              << (build_rotor ? rotor_keeper->getBaseParticles().TotalRealParticles() : 0) << " n_glass="
              << (build_glass ? glass_keeper->getBaseParticles().TotalRealParticles() : 0) << " n_wall="
              << (build_wall ? wall_keeper->getBaseParticles().TotalRealParticles() : 0)
              << " dp=" << stirring.french.dp << " steps=" << stirring.relax_steps << std::endl;

    for (int ite = 0; ite < stirring.relax_steps; ++ite)
    {
        update_configuration.exec();
        relaxation_residual.exec();
        const Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);

        if (stirring.relax_vtp_every > 0 && ite > 0 && (ite % stirring.relax_vtp_every == 0))
        {
            std::cout << "[ophelie][stirring-relax] step " << ite << std::endl;
            body_state_recorder.writeToFile(ite);
        }
    }

    body_state_recorder.writeToFile(stirring.relax_steps);
    write_particle_reload_files.writeToFile(0);
    // Reload.xml holds every body in one file and is rewritten from scratch, so a partial
    // --bodies run drops the body it did not build.
    std::cout << "[ophelie][stirring-relax] done; " << IO::getEnvironment().ReloadFolder()
              << "/Reload.xml now contains bodies=" << stirring.relax_bodies
              << (stirring.relax_bodies == "all" ? "" : " ONLY (rerun with --bodies=all before the flow case)")
              << std::endl;
    return 0;
#endif
}
