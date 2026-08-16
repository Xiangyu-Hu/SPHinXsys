/**
 * @file test_3d_ophelie_french_natural_convection_frozen_q.cpp
 * @brief Stage 4.1: freeze Q_50kW from EM+A_glass, then WCSPH + Boussinesq + Natural thermal BC.
 *
 * No EM update during flow. Side/bottom no-slip wall; open top ≈ free-slip.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

inline bool reloadXmlExists(const std::string &folder)
{
    return fs::exists(fs::path(folder) / "Reload.xml");
}

inline std::string resolveDefaultFrenchReloadFolder()
{
    const StdVec<std::string> candidates = {
        "./reload",
        "../reload",
        "../../../../../reload",
        "../test_3d_ophelie_french_natural_glass_relax/bin/reload",
        "../../test_3d_ophelie_french_natural_glass_relax/bin/reload",
    };
    for (const std::string &candidate : candidates)
    {
        if (reloadXmlExists(candidate))
        {
            return candidate;
        }
    }
    return candidates.front();
}

struct LocalCli
{
    std::string reload_dir;
    Real end_time = 0.02;
    Real u_ref = 0.05;
    Real c0 = 5.0;
    Real beta = 1.0e-5;
    Real gravity_g = 9.81;
    Real wall_thickness_factor = 2.0;
    bool enable_self_induction = true;
    bool enable_target_power = true;
    Real target_power_w = 50000.0;
};

inline void applyLocalCli(int ac, char *av[], OphelieFrenchReducedCaseParams &french, LocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--end-time=", 11) == 0)
        {
            cli.end_time = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--u-ref=", 8) == 0)
        {
            cli.u_ref = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--c0=", 5) == 0)
        {
            cli.c0 = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--beta=", 7) == 0)
        {
            cli.beta = static_cast<Real>(std::atof(av[i] + 7));
        }
        else if (std::strncmp(av[i], "--target-power=", 15) == 0)
        {
            cli.enable_target_power = true;
            cli.target_power_w = static_cast<Real>(std::atof(av[i] + 15));
        }
        else if (std::strcmp(av[i], "--no-self-induction") == 0)
        {
            cli.enable_self_induction = false;
        }
    }
    if (cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultFrenchReloadFolder();
    }
}

/** Explicit heating that keeps OphelieThermalDeltaT synchronized with Temperature. */
class ApplyFrozenQToTemperatureCK : public LocalDynamics
{
  public:
    ApplyFrozenQToTemperatureCK(SPHBody &sph_body, Real rho, Real cp, Real t_initial)
        : LocalDynamics(sph_body), inv_rho_cp_(Real(1) / (rho * cp + TinyReal)), t_initial_(t_initial),
          dv_q_(particles_->template getVariableByName<Real>("JouleHeat")),
          dv_temperature_(particles_->template getVariableByName<Real>(kOphelieTemperatureField)),
          dv_delta_t_(particles_->template getVariableByName<Real>(kOphelieThermalDeltaTField))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : inv_rho_cp_(encloser.inv_rho_cp_), t_initial_(encloser.t_initial_),
              q_(encloser.dv_q_->DelegatedData(ex_policy)),
              temperature_(encloser.dv_temperature_->DelegatedData(ex_policy)),
              delta_t_(encloser.dv_delta_t_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt)
        {
            temperature_[index_i] += q_[index_i] * dt * inv_rho_cp_;
            delta_t_[index_i] = temperature_[index_i] - t_initial_;
        }

      protected:
        Real inv_rho_cp_;
        Real t_initial_;
        Real *q_;
        Real *temperature_;
        Real *delta_t_;
    };

  protected:
    Real inv_rho_cp_;
    Real t_initial_;
    DiscreteVariable<Real> *dv_q_;
    DiscreteVariable<Real> *dv_temperature_;
    DiscreteVariable<Real> *dv_delta_t_;
};

#if SPHINXSYS_USE_SYCL
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
}
#else
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
}
#endif

} // namespace

int main(int ac, char *av[])
{
    LocalCli local_cli;
    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    // Natural geometry defaults (override reduced CAD placeholders).
    french.glass_radius = 0.25;
    french.glass_half_height = 0.5 * 0.185;
    french.glass_center = Vecd(0.0, 0.0, french.glass_half_height);
    french.dp = 0.015;
    french.coil.loop_radius = 0.285;
    french.coil.num_loops = 7;
    french.coil.segments_per_loop = 64;
    french.coil.z_min = -0.0225;
    french.coil.z_max = 0.2075;
    french.auto_coil_z = false;
    params.frequency_ = 282000.0;
    params.sigma_glass_ = 16.0;

    applyLocalCli(ac, av, french, local_cli);
    const StdVec<std::string> french_filtered = filterFrenchReducedCommandLine(ac, av, french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    OphelieTestCliOptions cli_options;
    StdVec<char *> french_av;
    for (auto &argument : french_filtered)
    {
        french_av.push_back(const_cast<char *>(argument.c_str()));
    }
    (void)filterOphelieTestCommandLine(static_cast<int>(french_av.size()), french_av.data(), params, cli_options);
    if (!cli_options.reload_dir.empty())
    {
        local_cli.reload_dir = cli_options.reload_dir;
    }

    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_self_induction_ = local_cli.enable_self_induction;
    if (local_cli.enable_target_power)
    {
        params.target_joule_power_ = local_cli.target_power_w;
        params.enable_power_scaling_ = false;
    }
    if (params.self_induction_phi_eq_res_tolerance_ > Real(1.0e-3))
    {
        params.self_induction_phi_eq_res_tolerance_ = Real(2.0e-4);
    }
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);
    logOphelieFinalParams(params, cli_options);

    if (!reloadXmlExists(local_cli.reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_natural_convection_frozen_q: Reload.xml required under \""
                  << local_cli.reload_dir << "\"\n";
        return 1;
    }

    //----------------------------------------------------------------------
    // Phase A — EM on SolidBody; freeze Q_50kW particle field.
    // Keep EM SPHSystem alive for process lifetime: tearing down a SYCL
    // SPHSystem then creating another often hangs the GPU runtime.
    //----------------------------------------------------------------------
    StdVec<Real> frozen_q;
    Real p_joule = 0.0;
    Real phi_eq = 0.0;
    UniquePtr<SPHSystem> em_system;
    UniquePtr<SolidBody> glass_em;
    UniquePtr<Inner<>> glass_em_inner;
    UniquePtr<RegisterOphelieGlassFields> register_glass_em;
    {
        const BoundingBoxd bounds = frenchReducedDomainBounds(french, 3.0 * french.dp);
        em_system = makeUnique<SPHSystem>(bounds, french.dp);
        em_system->setReloadParticles(true);
        IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

        glass_em = makeUnique<SolidBody>(
            *em_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                            "GlassBody", french.glass_center, french.glass_radius, french.glass_half_height));
        glass_em->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        glass_em->defineMatterMaterial<Solid>();
        (void)defineOphelieSolidLevelSet(*glass_em);
        glass_em->generateParticles<BaseParticles, Reload>(glass_em->Name());
        em_system->initializeSystemCellLinkedLists();
        em_system->initializeSystemConfigurations();

        OphelieGlassFieldNames glass_names;
        register_glass_em = makeUnique<RegisterOphelieGlassFields>(*glass_em, glass_names);
        glass_em_inner = makeUnique<Inner<>>(*glass_em);
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(*glass_em, glass_names,
                                                                                    params.sigma_glass_);
        assign_sigma.exec();

        OphelieFrenchEmJouleHeatOneWayResult probe;
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            *glass_em, *glass_em_inner, glass_names, params, french, probe);
        if (local_cli.enable_target_power)
        {
            (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
            runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                *glass_em, *glass_em_inner, glass_names, params, french, probe);
        }
        p_joule = probe.joule_power_w;
        phi_eq = probe.phi_eq_res_vol;

        BaseParticles &particles = glass_em->getBaseParticles();
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_names, params);
        syncVariableToHost<Real>(particles, q_field);
        const Real *q = particles.getVariableDataByName<Real>(q_field);
        const size_t n = particles.TotalRealParticles();
        frozen_q.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            frozen_q[i] = q[i];
        }
        std::cout << "[ophelie][stage4.1] EM freeze: n=" << n << " P_joule_W=" << p_joule
                  << " phi_eq_res_vol=" << phi_eq << " self_induction=" << (probe.used_self_induction ? 1 : 0)
                  << std::endl;
    }

    //----------------------------------------------------------------------
    // Phase B — FluidBody WCSPH + Boussinesq + Natural thermal BC.
    // Domain is glass±wall only (not coil): Lattice scans the full SPHSystem
    // AABB, so a coil-sized box makes wall generation extremely slow.
    //----------------------------------------------------------------------
    const Real wall_thickness = local_cli.wall_thickness_factor * french.dp;
    const Real flow_margin = Real(4.0) * french.dp + wall_thickness;
    const Vecd &gc = french.glass_center;
    const Real flow_r = french.glass_radius + wall_thickness + flow_margin;
    const Real flow_hz = french.glass_half_height + wall_thickness + flow_margin;
    const BoundingBoxd flow_bounds(Vecd(gc[0] - flow_r, gc[1] - flow_r, gc[2] - flow_hz),
                                   Vecd(gc[0] + flow_r, gc[1] + flow_r, gc[2] + flow_hz));
    SPHSystem sph_system(flow_bounds, french.dp);
    sph_system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

    const OphelieFrenchGlassMaterialLaws laws = makeFrenchNaturalCep2008MaterialLaws();
    const Real t0 = Real(1473);
    const Real rho0 = laws.rho.evaluate(t0);
    const Real mu = laws.mu.evaluate(t0);
    const Real cp = laws.cp.evaluate(t0);
    const Real k_th = laws.k.evaluate(t0);
    const Real beta = local_cli.beta > TinyReal ? local_cli.beta : laws.beta.evaluate(t0);

    FluidBody glass(sph_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                                    "GlassBody", french.glass_center, french.glass_radius, french.glass_half_height));
    glass.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass.defineMatterMaterial<WeaklyCompressibleFluid>(rho0, local_cli.c0);
    glass.addMaterialProperty<Viscosity>(mu);
    glass.generateParticles<BaseParticles, Reload>(glass.Name());

    // SPHinXsys cylinder wall path: TriangleMeshShapeCylinder + LevelSet + NormalFromBodyShapeCK
    // (same pattern as taylor_bar_sycl / french_natural_glass_relax).
    const int wall_mesh_resolution =
        french.glass_mesh_resolution > 0 ? french.glass_mesh_resolution : 20;
    SolidBody wall(sph_system, makeShared<OphelieFrenchNaturalCrucibleWallShape>(
                                   "WallBoundary", french, wall_thickness, wall_mesh_resolution));
    wall.defineAdaptationRatios(1.3, 1.0);
    wall.defineMatterMaterial<Solid>();
    (void)wall.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
    wall.generateParticles<BaseParticles, Lattice>();

    Inner<> glass_inner(glass);
    Contact<> glass_wall_contact(glass, {&wall});

    BaseParticles &glass_particles = glass.getBaseParticles();
    registerOphelieJouleHeatTemperatureField(glass_particles, t0);
    registerOphelieThermalDiffusionAuxFields(glass_particles, k_th);
    glass_particles.registerStateVariable<Real>("JouleHeat", Real(0));
    {
        Real *q_host = glass_particles.getVariableDataByName<Real>("JouleHeat");
        const size_t n = glass_particles.TotalRealParticles();
        if (n != frozen_q.size())
        {
            std::cerr << "[ophelie][stage4.1] particle count mismatch EM n=" << frozen_q.size() << " fluid n=" << n
                      << std::endl;
            return 1;
        }
        for (size_t i = 0; i < n; ++i)
        {
            q_host[i] = frozen_q[i];
        }
        syncVariableToDevice<Real>(glass_particles, "JouleHeat");
    }

    OphelieThermalDiffusionOneWayOptions thermal_bc;
    thermal_bc.enable_french_natural_bc = true;
    thermal_bc.enable_cold_wall_dirichlet = false;
    thermal_bc.enable_diffusion = false;
    thermal_bc.boundary_width_factor = params.phi_boundary_distance_factor_;
    thermal_bc.h_side = Real(300);
    thermal_bc.h_bottom = Real(35);
    thermal_bc.h_free = Real(20);
    thermal_bc.emissivity = Real(0.8);
    thermal_bc.t_cool = Real(300);
    thermal_bc.t_ambient = Real(300);
    thermal_bc.t_rad_ambient = Real(300);
    (void)setupOphelieThermalFrenchNaturalBoundaryFaces(glass_particles, thermal_bc, french);

    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall).exec();

    auto &update_glass_cell_linked_list = main_methods.addCellLinkedListDynamics(glass);
    auto &update_wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall);
    auto &update_glass_relations = main_methods.addRelationDynamics(glass_inner, glass_wall_contact);
    auto &fluid_particle_sorting = main_methods.addSortDynamics(glass);

    Gravity gravity(Vecd(0.0, 0.0, -local_cli.gravity_g));
    auto &constant_gravity = main_methods.addStateDynamics<GravityForceCK<Gravity>>(glass, gravity);
    auto &boussinesq_force =
        main_methods.addStateDynamics<BoussinesqBuoyancyForceCK>(glass, gravity, beta, t0, kOphelieTemperatureField);
    auto &fluid_advection_step_setup = main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(glass);
    auto &fluid_update_particle_position = main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(glass);

    auto &acoustic_step_1st_half =
        main_methods
            .addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep1stHalf, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>(glass_wall_contact);
    auto &acoustic_step_2nd_half =
        main_methods.addInteractionDynamicsOneLevel<fluid_dynamics::AcousticStep2ndHalf, AcousticRiemannSolverCK,
                                                    NoKernelCorrectionCK>(glass_inner);
    auto &acoustic_step_2nd_half_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::AcousticStep2ndHalf, Wall, AcousticRiemannSolverCK,
                                            NoKernelCorrectionCK>(glass_wall_contact);
    acoustic_step_2nd_half.addPostContactInteraction(acoustic_step_2nd_half_wall);
    auto &density_regularization =
        main_methods.addInteractionDynamics<fluid_dynamics::DensitySummationCK>(glass_inner)
            .addPostContactInteraction(glass_wall_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(glass);

    auto &fluid_viscous_force =
        main_methods.addInteractionDynamicsWithUpdate<fluid_dynamics::ViscousForceCK, Viscosity, NoKernelCorrectionCK>(
            glass_inner);
    auto &fluid_viscous_force_from_wall =
        main_methods.addInteractionDynamics<fluid_dynamics::ViscousForceCK, Wall, Viscosity, NoKernelCorrectionCK>(
            glass_wall_contact);
    fluid_viscous_force.addPostContactInteraction(fluid_viscous_force_from_wall);

    auto &advection_time_step = main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(glass, local_cli.u_ref);
    auto &acoustic_time_step = main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(glass, 0.55);
    auto &u_max_reduce = main_methods.addReduceDynamics<OphelieVelocityMaxReduceCK>(glass);
    auto &t_mean_reduce = main_methods.addReduceDynamics<OphelieTemperatureMeanReduceCK>(glass);
    auto &t_max_reduce =
        main_methods.addReduceDynamics<OphelieThermalMaxTemperatureReduceCK>(glass, kOphelieTemperatureField);

    StateDynamics<MainExecutionPolicy, ApplyFrozenQToTemperatureCK> apply_q(glass, rho0, cp, t0);
    StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_natural_bc(
        glass, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, thermal_bc.boundary_shell_thickness,
        thermal_bc);

    update_wall_cell_linked_list.exec();
    update_glass_cell_linked_list.exec();
    update_glass_relations.exec();
    density_regularization.exec();
    fluid_advection_step_setup.exec();
    constant_gravity.exec();
    boussinesq_force.exec();

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    size_t advection_steps = 0;

    std::cout << "[ophelie][stage4.1] flow start: end_time=" << local_cli.end_time
              << " n_glass=" << glass_particles.TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " beta=" << beta << " mu=" << mu
              << " c0=" << local_cli.c0 << " T0=" << t0 << std::endl;

    while (!time_stepper.isEndTime(local_cli.end_time))
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);

        apply_q.exec(acoustic_dt);
        apply_natural_bc.exec(acoustic_dt);
        boussinesq_force.exec();

        if (advection_step(advection_time_step))
        {
            ++advection_steps;
            fluid_update_particle_position.exec();
            if (advection_steps % 50 == 0)
            {
                fluid_particle_sorting.exec();
            }
            update_glass_cell_linked_list.exec();
            update_glass_relations.exec();
            density_regularization.exec();
            fluid_advection_step_setup.exec();
            fluid_viscous_force.exec();
            constant_gravity.exec();
            boussinesq_force.exec();
            if (advection_steps % 20 == 0)
            {
                std::cout << "[ophelie][stage4.1] N=" << advection_steps << " t=" << time_stepper.getPhysicalTime()
                          << " U_max=" << u_max_reduce.exec() << " T_mean=" << t_mean_reduce.exec()
                          << " T_max=" << t_max_reduce.exec() << std::endl;
            }
        }
    }

    const Real u_max = u_max_reduce.exec();
    const Real t_mean = t_mean_reduce.exec();
    const Real t_max = t_max_reduce.exec();
    Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                  side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);

    const Real power_rel_err =
        std::abs(p_joule - local_cli.target_power_w) / (local_cli.target_power_w + TinyReal);
    const bool em_ok = std::isfinite(p_joule) && power_rel_err < Real(1.0e-2);
    const bool flow_ok = std::isfinite(u_max) && std::isfinite(t_mean) && std::isfinite(t_max) && u_max >= Real(0);
    const bool thermal_ok = total_loss_w > TinyReal && side_w > TinyReal && bottom_w > TinyReal &&
                            (free_conv_w + free_rad_w) > TinyReal;
    const bool buoyancy_ok = u_max > TinyReal; // expect some motion from Boussinesq + heating/cooling
    const bool passed = em_ok && flow_ok && thermal_ok && buoyancy_ok;

    std::cout << "test_3d_ophelie_french_natural_convection_frozen_q"
              << " P_joule_W=" << p_joule << " power_rel_err=" << power_rel_err << " phi_eq_res_vol=" << phi_eq
              << " end_time=" << local_cli.end_time << " U_max=" << u_max << " T_mean=" << t_mean << " T_max=" << t_max
              << " wall_loss_side_W=" << side_w << " wall_loss_bottom_W=" << bottom_w
              << " free_conv_loss_W=" << free_conv_w << " free_rad_loss_W=" << free_rad_w
              << " total_heat_loss_W=" << total_loss_w << " em_ok=" << (em_ok ? 1 : 0)
              << " flow_ok=" << (flow_ok ? 1 : 0) << " thermal_ok=" << (thermal_ok ? 1 : 0)
              << " buoyancy_ok=" << (buoyancy_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
