/**
 * @file test_3d_ophelie_client_nc_glass_relax.cpp
 * @brief Client natural-convection parametric melt: lattice + SYCL-CK particle relaxation.
 *
 * Geometry follows Doc software plan / Doc229:
 *   - melt from inner radius + provisional liquid depth (NOT shell STL fill)
 *   - R≈0.1415 m, H_melt default 0.180 m (provisional), bottom at z≈0.001 m
 *   - writes melt/probe JSON and effective material snapshot; Reload for later NC short runs
 *
 * Pattern mirrors test_3d_ophelie_french_natural_glass_relax.
 */
#include "electromagnetic_ophelie_client_case_params.h"
#include "electromagnetic_ophelie_client_geometry.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
namespace fs = std::filesystem;

namespace
{
OphelieClientMeltGeometryM g_melt = makeClientNcMeltGeometryDevelopment();
Real particle_spacing_ref = Real(0.010);
int resolution = 20;
int relaxation_steps = 1000;
int relax_vtp_every = 100;
std::string output_meta_dir = "client_nc_relax_meta";
bool use_client_ln_sigma = false;
BoundingBoxd system_domain_bounds(Vecd::Zero(), Vecd::Zero());

inline void refreshDomain()
{
    const Real pad = Real(4) * particle_spacing_ref;
    system_domain_bounds = clientNcDomainBounds(g_melt, pad);
}

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
        if (std::strncmp(arg, "--melt-height=", 14) == 0)
        {
            g_melt.melt_height_m = static_cast<Real>(std::atof(arg + 14));
            g_melt.glass_half_height = Real(0.5) * g_melt.melt_height_m;
            g_melt.glass_center_z_m = g_melt.inner_bottom_z_m + g_melt.glass_half_height;
            g_melt.glass_center = Vecd(0, 0, g_melt.glass_center_z_m);
            g_melt.melt_height_status = "user_cli_override";
            adaptClientNcCoilStackToMeltHeight(g_melt);
            continue;
        }
        if (std::strncmp(arg, "--melt-radius=", 14) == 0)
        {
            g_melt.melt_radius_m = static_cast<Real>(std::atof(arg + 14));
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
        if (std::strncmp(arg, "--output-meta-dir=", 18) == 0)
        {
            output_meta_dir = arg + 18;
            continue;
        }
        if (std::strcmp(arg, "--client-ln-sigma") == 0)
        {
            use_client_ln_sigma = true;
            continue;
        }
        sph_av.push_back(const_cast<char *>(arg));
    }
    refreshDomain();
}

inline void writeMetaAndMaterial()
{
    fs::create_directories(output_meta_dir);
    const auto cat = makeClientNaturalConvectionAssemblyCatalog();
    const auto probe = makeClientNominalThermocoupleProbe(cat);
    writeClientNcMeltGeometryJson(g_melt, probe, (fs::path(output_meta_dir) / "melt_geometry.json").string());

    OphelieClientCaseParams preset = use_client_ln_sigma ? makeClientNaturalConvectionClientSigmaContrastParams()
                                                         : makeClientNaturalConvectionDevelopmentParams();
    preset.melt_height_m = g_melt.melt_height_m;
    preset.particle_spacing_m = particle_spacing_ref;
    writeClientEffectiveCaseJson(preset, (fs::path(output_meta_dir) / "effective-case.json").string());

    const auto laws = buildClientCaseMaterialLaws(preset);
    logClientMaterialIdentity(preset.material_id, laws, preset.mu_freeze_t_k);

    std::ofstream probe_csv(fs::path(output_meta_dir) / "probe_nominal.csv");
    probe_csv << "probe_id,status,x_m,y_m,z_m,target_T_K,note\n";
    probe_csv << std::setprecision(10) << probe.probe_id << ',' << probe.status << ',' << probe.solver_m.x << ','
              << probe.solver_m.y << ',' << probe.solver_m.z << ',' << probe.target_temperature_k << ','
              << "\"nominal AABB tip; unconfirmed sensor\"\n";

    std::ofstream readme(fs::path(output_meta_dir) / "README.md");
    readme << "# Client NC glass relax meta\n\n";
    readme << "- parametric melt R=" << g_melt.melt_radius_m << " H=" << g_melt.melt_height_m << " ("
           << g_melt.melt_height_status << ")\n";
    readme << "- dp=" << particle_spacing_ref << "\n";
    readme << "- material_id=" << preset.material_id << "\n";
    readme << "- frozen mu @ " << preset.mu_freeze_t_k << " K = " << laws.mu.evaluate(preset.mu_freeze_t_k) << "\n";
    readme << "- sigma @ ref = " << laws.sigma.evaluate(preset.mu_freeze_t_k) << "\n";
    readme << "- generator_input_power_W=30000 (not absorbed)\n";
    readme << "- NOT a 1000C quantitative PASS\n";
}
} // namespace

int main(int ac, char *av[])
{
    StdVec<char *> sph_av;
    applyLocalCommandLine(ac, av, sph_av);

    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(static_cast<int>(sph_av.size()), sph_av.data());

    std::cout << "[ophelie][client-nc-relax] R=" << g_melt.melt_radius_m << " H=" << g_melt.melt_height_m
              << " (" << g_melt.melt_height_status << ") center_z=" << g_melt.glass_center[2]
              << " dp=" << particle_spacing_ref << " mesh_resolution=" << resolution
              << " relax_steps=" << relaxation_steps << std::endl;

    writeMetaAndMaterial();

    RealBody glass_body(sph_system, makeShared<OphelieClientNcGlassCylinderShape>("GlassBody", g_melt, resolution));
    LevelSetShape &level_set_shape = glass_body.defineBodyLevelSetShape(par_ck, 2.0).correctLevelSetSign().writeLevelSet();
    glass_body.defineMatterMaterial<Solid>();
    glass_body.generateParticles<BaseParticles, Lattice>();

    NearShapeSurface near_body_surface(glass_body);
    Inner<> glass_body_inner(glass_body);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    auto &glass_body_cell_linked_list = main_methods.addCellLinkedListDynamics(glass_body);
    auto &glass_body_update_inner_relation = main_methods.addRelationDynamics(glass_body_inner);
    auto &random_glass_body_particles = host_methods.addStateDynamics<RandomizeParticlePositionCK>(glass_body);
    auto &relaxation_residual =
        main_methods.addInteractionDynamics<KernelGradientIntegral, NoKernelCorrectionCK>(glass_body_inner)
            .addPostStateDynamics<LevelsetKernelGradientIntegral>(glass_body, level_set_shape);
    auto &relaxation_scaling = main_methods.addReduceDynamics<RelaxationScalingCK>(glass_body);
    auto &update_particle_position = main_methods.addStateDynamics<PositionRelaxationCK>(glass_body);
    auto &level_set_bounding = main_methods.addStateDynamics<LevelsetBounding>(near_body_surface);

    auto &body_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass_body);
    auto &write_particle_reload_files =
        main_methods.addIODynamics<ReloadParticleIOCK>(StdVec<SPHBody *>{&glass_body});

    random_glass_body_particles.exec();
    glass_body_cell_linked_list.exec();
    level_set_bounding.exec();
    if (relax_vtp_every > 0)
        body_state_recorder.writeToFile(0);

    int ite_p = 0;
    while (ite_p < relaxation_steps)
    {
        glass_body_cell_linked_list.exec();
        glass_body_update_inner_relation.exec();
        relaxation_residual.exec();
        Real relaxation_step = relaxation_scaling.exec();
        update_particle_position.exec(relaxation_step);
        level_set_bounding.exec();
        ++ite_p;
        if (relax_vtp_every > 0 && ite_p % relax_vtp_every == 0)
        {
            body_state_recorder.writeToFile(ite_p);
            std::cout << "[ophelie][client-nc-relax] step " << ite_p << " / " << relaxation_steps << std::endl;
        }
    }

    write_particle_reload_files.writeToFile(0);
    std::cout << "[ophelie][client-nc-relax] done. meta=" << output_meta_dir
              << " N_real=" << glass_body.getBaseParticles().TotalRealParticles() << std::endl;
    return 0;
}
