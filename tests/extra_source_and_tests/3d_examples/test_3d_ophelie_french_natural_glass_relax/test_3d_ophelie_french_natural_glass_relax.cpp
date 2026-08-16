/**
 * @file test_3d_ophelie_french_natural_glass_relax.cpp
 * @brief French natural-convection glass cylinder: generate geometry, SYCL level-set,
 *        lattice particles, then SYCL-CK particle relaxation.
 *
 * Geometry defaults (Jacoutot natural case / EREBUS fill):
 *   R = 0.250 m, H = 0.185 m, axis +z, z in [0, H].
 *
 * Level-set and relaxation follow the official SPHinXsys SYCL pattern in
 * tests/tests_sycl/3d_examples/test_3d_particle_relaxation_single_resolution_sycl
 * and test_3d_taylor_bar_sycl.
 */
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>

using namespace SPH;

//----------------------------------------------------------------------
//	Runtime / default geometry parameters.
//----------------------------------------------------------------------
Real glass_radius = 0.250;
Real glass_height = 0.185;
// Stage1.5: align production dp with reload spacing used for French natural.
Real particle_spacing_ref = 0.015;
int resolution = 20;
int relaxation_steps = 1000;
int relax_vtp_every = 100;
Vec3d translation_column(0.0, 0.0, 0.5 * glass_height);
BoundingBoxd system_domain_bounds(Vec3d::Zero(), Vec3d::Zero());

inline void refreshGeometryDependentQuantities()
{
    translation_column = Vec3d(0.0, 0.0, 0.5 * glass_height);
    const Real pad = 4.0 * particle_spacing_ref;
    system_domain_bounds =
        BoundingBoxd(Vec3d(-glass_radius - pad, -glass_radius - pad, -pad),
                     Vec3d(glass_radius + pad, glass_radius + pad, glass_height + pad));
}

/** French natural glass cylinder. */
class GlassBody : public ComplexShape
{
  public:
    explicit GlassBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeCylinder>(Vec3d(0.0, 0.0, 1.0), glass_radius, 0.5 * glass_height, resolution,
                                       translation_column);
    }
};

inline void applyLocalCommandLine(int ac, char *av[], StdVec<char *> &sph_av)
{
    sph_av.push_back(av[0]);
    for (int i = 1; i < ac; ++i)
    {
        const char *arg = av[i];
        if (std::strncmp(arg, "--dp=", 5) == 0)
        {
            particle_spacing_ref = static_cast<Real>(std::atof(arg + 5));
            continue;
        }
        if (std::strncmp(arg, "--glass-radius=", 15) == 0)
        {
            glass_radius = static_cast<Real>(std::atof(arg + 15));
            continue;
        }
        if (std::strncmp(arg, "--glass-height=", 15) == 0)
        {
            glass_height = static_cast<Real>(std::atof(arg + 15));
            continue;
        }
        if (std::strncmp(arg, "--mesh-resolution=", 18) == 0)
        {
            resolution = std::atoi(arg + 18);
            continue;
        }
        if (std::strncmp(arg, "--relax-steps=", 14) == 0)
        {
            relaxation_steps = std::atoi(arg + 14);
            continue;
        }
        if (std::strncmp(arg, "--relax-vtp-every=", 18) == 0)
        {
            relax_vtp_every = std::atoi(arg + 18);
            continue;
        }
        sph_av.push_back(const_cast<char *>(arg));
    }
    refreshGeometryDependentQuantities();
}

int main(int ac, char *av[])
{
    StdVec<char *> sph_av;
    applyLocalCommandLine(ac, av, sph_av);
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());

    std::cout << "French natural glass relax: R=" << glass_radius << " H=" << glass_height
              << " dp=" << particle_spacing_ref << " mesh_resolution=" << resolution
              << " relax_steps=" << relaxation_steps << " vtp_every=" << relax_vtp_every << std::endl;
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    RealBody glass_body(sph_system, makeShared<GlassBody>("GlassBody"));
    // level set shape is used for particle relaxation (SYCL)
    // refinement=2.0 matches tests/tests_sycl/3d_examples/test_3d_taylor_bar_sycl
    LevelSetShape &level_set_shape = glass_body.defineBodyLevelSetShape(par_ck, 2.0)
                                         .correctLevelSetSign()
                                         .writeLevelSet();
    glass_body.defineMatterMaterial<Solid>();
    glass_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Creating body parts.
    //----------------------------------------------------------------------
    NearShapeSurface near_body_surface(glass_body);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    Inner<> glass_body_inner(glass_body);
    //----------------------------------------------------------------------
    //	Methods used for particle relaxation.
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    //----------------------------------------------------------------------
    auto &glass_body_cell_linked_list = main_methods.addCellLinkedListDynamics(glass_body);
    auto &glass_body_update_inner_relation = main_methods.addRelationDynamics(glass_body_inner);
    auto &random_glass_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(glass_body);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(glass_body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(glass_body, level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(glass_body);
    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(glass_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass_body);
    auto &write_particle_reload_files =
        main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&glass_body});
    //----------------------------------------------------------------------
    //	Prepare: randomize, then surface bounding before the first VTP.
    //	Official CPU taylor_bar does SurfaceBounding before writeToFile(0).
    //----------------------------------------------------------------------
    random_glass_body_particles.exec();
    glass_body_cell_linked_list.exec();
    level_set_bounding.exec();
    //----------------------------------------------------------------------
    //	First output before the simulation.
    //----------------------------------------------------------------------
    if (relax_vtp_every > 0)
    {
        body_state_recorder.writeToFile(0);
    }
    //----------------------------------------------------------------------
    //	Particle relaxation time stepping start here.
    //----------------------------------------------------------------------
    int ite_p = 0;
    while (ite_p < relaxation_steps)
    {
        glass_body_cell_linked_list.exec();
        glass_body_update_inner_relation.exec();

        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();

        ite_p += 1;
        if (ite_p % 100 == 0 || ite_p == relaxation_steps)
        {
            std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
        }
        if (relax_vtp_every > 0 && (ite_p % relax_vtp_every == 0 || ite_p == relaxation_steps))
        {
            body_state_recorder.writeToFile(ite_p);
        }
    }
    std::cout << "The physics relaxation process finish !" << std::endl;
    write_particle_reload_files.writeToFile();

    return 0;
}
