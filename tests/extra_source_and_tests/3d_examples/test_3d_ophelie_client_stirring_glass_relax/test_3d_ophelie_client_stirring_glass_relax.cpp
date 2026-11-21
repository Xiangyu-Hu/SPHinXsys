/**
 * @file test_3d_ophelie_client_stirring_glass_relax.cpp
 * @brief Client stirring relax — same body layout as French stirring
 *        (GlassBody / Rotor / WallBoundary → Reload), client geometry + materials.
 *
 * Glass = parametric melt cylinder minus paddle STL (shell STL display-only).
 * Rotor = preprocessed solver-frame paddle (metres, +Z up).
 * Wall  = numerical open cup (OphelieClientNcCrucibleWallShape).
 *
 * Next: elastic_driven Rotor in a french_stirring_em-structured flow/EM case.
 */
#include "electromagnetic_ophelie_client_geometry.h"
#include "electromagnetic_ophelie_client_stirring_geometry.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cstring>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;

namespace
{
OphelieClientStirringCaseParams g_case;
bool build_rotor = true;
bool build_glass = true;
bool build_wall = true;

void applyCli(int ac, char *av[], StdVec<char *> &sph_av)
{
    sph_av.push_back(av[0]);
    for (int i = 1; i < ac; ++i)
    {
        const char *arg = av[i];
        if (std::strncmp(arg, "--rotor-stl=", 12) == 0)
        {
            g_case.rotor_stl_solver_m = arg + 12;
            continue;
        }
        if (std::strncmp(arg, "--dp=", 5) == 0)
        {
            g_case.client.particle_spacing_m = static_cast<Real>(std::atof(arg + 5));
            continue;
        }
        if (std::strncmp(arg, "--melt-height=", 14) == 0)
        {
            g_case.melt_height_m = static_cast<Real>(std::atof(arg + 14));
            continue;
        }
        if (std::strncmp(arg, "--relax-steps=", 14) == 0)
        {
            g_case.french_stirring.relax_steps = std::atoi(arg + 14);
            continue;
        }
        if (std::strncmp(arg, "--relax-vtp-every=", 18) == 0)
        {
            g_case.french_stirring.relax_vtp_every = std::atoi(arg + 18);
            continue;
        }
        if (std::strncmp(arg, "--bodies=", 9) == 0)
        {
            g_case.french_stirring.relax_bodies = arg + 9;
            const std::string s = g_case.french_stirring.relax_bodies;
            if (s == "all")
            {
                build_rotor = build_glass = build_wall = true;
            }
            else
            {
                build_rotor = s.find("rotor") != std::string::npos;
                build_glass = s.find("glass") != std::string::npos;
                build_wall = s.find("wall") != std::string::npos;
            }
            continue;
        }
        sph_av.push_back(const_cast<char *>(arg));
    }
    applyClientStirringDevelopmentDefaults(g_case);
}
} // namespace

int main(int ac, char *av[])
{
#if !SPHINXSYS_USE_SYCL
    (void)ac;
    (void)av;
    std::cerr << "test_3d_ophelie_client_stirring_glass_relax requires SYCL build.\n";
    return 1;
#else
    StdVec<char *> sph_av;
    applyCli(ac, av, sph_av);
    applyClientStirringDevelopmentDefaults(g_case);

    if (!build_rotor && !build_glass && !build_wall)
    {
        std::cerr << "[ophelie][client-stir-relax] --bodies must be all or subset of glass,rotor,wall\n";
        return 2;
    }

    const BoundingBoxd system_domain_bounds = clientStirringRelaxDomainBounds(g_case);
    SPHSystem sph_system(system_domain_bounds, g_case.french_stirring.french.dp);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());
    IO::getEnvironment().resetReloadFolder("./reload", true);

    std::cout << "[ophelie][client-stir-relax] mechanics=" << clientPaddleMechanicsName(g_case.paddle_mechanics)
              << " rotor_stl=" << g_case.rotor_stl_solver_m << " R=" << g_case.melt_radius_m
              << " H=" << g_case.melt_height_m << " dp=" << g_case.french_stirring.french.dp
              << " bodies=" << g_case.french_stirring.relax_bodies
              << " rpm=" << g_case.client.stirrer_rpm << "\n";

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
        rotor_keeper =
            makeUnique<SolidBody>(sph_system, makeShared<OphelieClientStirringRotorShape>("Rotor", g_case));
        SolidBody &rotor = *rotor_keeper;
        rotor.defineAdaptationRatios(1.3, 1.0);
        rotor.defineMatterMaterial<Solid>(g_case.steel.rho_kg_m3);
        LevelSetShape &rotor_ls =
            rotor.defineBodyLevelSetShape(par_ck, g_case.french_stirring.rotor_level_set_refinement)
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
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(rotor, rotor_ls));
        surface_bounding.add(&main_methods.addStateDynamics<LevelsetBounding>(*rotor_near_surface_keeper));
    }

    if (build_glass)
    {
        glass_keeper =
            makeUnique<SolidBody>(sph_system, makeShared<OphelieClientStirringGlassFluidShape>("GlassBody", g_case));
        SolidBody &glass = *glass_keeper;
        glass.defineAdaptation<SPHAdaptation>(1.0, 1.0);
        glass.defineMatterMaterial<Solid>();
        LevelSetShape &glass_ls =
            glass.defineBodyLevelSetShape(par_ck, g_case.french_stirring.glass_level_set_refinement)
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
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(glass, glass_ls));
        surface_bounding.add(&main_methods.addStateDynamics<LevelsetBounding>(*glass_near_surface_keeper));
    }

    if (build_wall)
    {
        OphelieClientMeltGeometryM melt = makeClientNcMeltGeometryDevelopment();
        melt.melt_height_m = g_case.melt_height_m;
        melt.glass_half_height = Real(0.5) * melt.melt_height_m;
        melt.glass_center_z_m = melt.inner_bottom_z_m + melt.glass_half_height;
        melt.glass_center = Vecd(0, 0, melt.glass_center_z_m);
        wall_keeper = makeUnique<SolidBody>(
            sph_system, makeShared<OphelieClientNcCrucibleWallShape>(
                            "WallBoundary", melt, Real(4) * g_case.french_stirring.french.dp, 20));
        SolidBody &wall = *wall_keeper;
        wall.defineAdaptationRatios(1.3, 1.0);
        wall.defineMatterMaterial<Solid>();
        LevelSetShape &wall_ls = wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
        wall.generateParticles<BaseParticles, Lattice>();
        wall_inner_keeper = makeUnique<Inner<>>(wall);
        wall_near_surface_keeper = makeUnique<NearShapeSurface>(wall);
        relax_bodies.push_back(&wall);
        reload_bodies.push_back(&wall);
        host_methods.addStateDynamics<RandomizeParticlePositionCK>(wall).exec();
        update_relation.add(&main_methods.addRelationDynamics(*wall_inner_keeper));
        relaxation_residual.add(
            &main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(*wall_inner_keeper)
                 .addPostStateDynamics<LevelsetKernelGradientIntegral>(wall, wall_ls));
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

    std::cout << "[ophelie][client-stir-relax] n_rotor="
              << (build_rotor ? rotor_keeper->getBaseParticles().TotalRealParticles() : 0) << " n_glass="
              << (build_glass ? glass_keeper->getBaseParticles().TotalRealParticles() : 0) << " n_wall="
              << (build_wall ? wall_keeper->getBaseParticles().TotalRealParticles() : 0)
              << " steps=" << g_case.french_stirring.relax_steps << "\n";

    for (int ite = 0; ite < g_case.french_stirring.relax_steps; ++ite)
    {
        update_configuration.exec();
        relaxation_residual.exec();
        const Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        if (g_case.french_stirring.relax_vtp_every > 0 && ite > 0 &&
            (ite % g_case.french_stirring.relax_vtp_every == 0))
        {
            std::cout << "[ophelie][client-stir-relax] step " << ite << "\n";
            body_state_recorder.writeToFile(ite);
        }
    }

    body_state_recorder.writeToFile(g_case.french_stirring.relax_steps);
    write_particle_reload_files.writeToFile(0);
    std::cout << "[ophelie][client-stir-relax] done; " << IO::getEnvironment().ReloadFolder()
              << "/Reload.xml bodies=" << g_case.french_stirring.relax_bodies << "\n";
    return 0;
#endif
}
