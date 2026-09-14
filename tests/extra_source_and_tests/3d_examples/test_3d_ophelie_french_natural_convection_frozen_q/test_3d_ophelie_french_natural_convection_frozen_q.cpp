/**
 * @file test_3d_ophelie_french_natural_convection_frozen_q.cpp
 * @brief Stage 4.1: frozen Q_50kW (default) or periodic thermo–EM coupling.
 *
 * Side/bottom no-slip wall with a short open rim above the melt (not a lid).
 * Hydrostatic density is applied, then a gravity-only settle, before Q / Boussinesq.
 * Periodic coupling defaults SPH conduction on (`--no-thermal-diffusion` to disable).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_boussinesq.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "electromagnetic_ophelie_thermo_em_coupling.h"
#include "io_environment.h"
#include "interaction_ck.h"
#include "sphinxsys.h"

#include "../test_3d_ophelie_rh200_glass_em_stirring/rh200_joule_heat_grid.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
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
    Real u_ref = 0.0;              // 0 -> sqrt(gH); was 0.05 and delayed the first advection
    Real c0 = 0.0;                 // 0 -> 10*sqrt(gH); was 5 and too soft for hydrostatics
    Real beta = 1.0e-5;
    Real gravity_g = 9.81;
    Real wall_thickness_factor = 4.0; // still-water / stirring dummy-wall width; 2 dp leaks through the floor
    Real wall_rim_height = -1.0;   // <0 -> 6*dp splash guard above the melt (still open top)
    Real hydro_settle_time = 0.15; // gravity-only settle before Q / Boussinesq
    Real escaped_frac_max = 0.01;
    Real state_record_interval = 0.0; // 0 -> no VTP
    bool thermal_diffusion_explicit = false;
    bool enable_self_induction = true;
    bool enable_target_power = true;
    Real target_power_w = 50000.0;
    FrenchLiteratureParameters literature = makeFrenchNaturalLiteratureParameters();
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
            cli.literature.enable_aind = false;
        }
        else if (tryApplyFrenchLiteratureCli(av[i], cli.literature, cli.thermal_diffusion_explicit))
        {
            if (std::strncmp(av[i], "--target-absorbed-power=", 24) == 0)
            {
                cli.enable_target_power = true;
                cli.target_power_w = cli.literature.target_glass_absorbed_power_w;
            }
            else if (std::strncmp(av[i], "--aind=", 7) == 0)
            {
                cli.enable_self_induction = cli.literature.enable_aind;
            }
            else if (std::strcmp(av[i], "--no-boussinesq") == 0)
            {
                // literature.enable_boussinesq already cleared
            }
        }
        else if (std::strncmp(av[i], "--wall-thickness-factor=", 24) == 0)
        {
            cli.wall_thickness_factor = static_cast<Real>(std::atof(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--wall-rim=", 11) == 0)
        {
            cli.wall_rim_height = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--hydro-settle-time=", 20) == 0)
        {
            cli.hydro_settle_time = static_cast<Real>(std::atof(av[i] + 20));
        }
        else if (std::strncmp(av[i], "--escaped-frac-max=", 19) == 0)
        {
            cli.escaped_frac_max = static_cast<Real>(std::atof(av[i] + 19));
        }
        else if (std::strncmp(av[i], "--state-record-interval=", 24) == 0)
        {
            cli.state_record_interval = static_cast<Real>(std::atof(av[i] + 24));
        }
    }
    if (cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultFrenchReloadFolder();
    }
    if (cli.enable_target_power)
    {
        cli.literature.target_glass_absorbed_power_w = cli.target_power_w;
    }
    finalizeFrenchPeriodicLiteratureDefaults(cli.literature, cli.thermal_diffusion_explicit);
}

inline void finalizeNaturalHydroDefaults(LocalCli &cli, const OphelieFrenchReducedCaseParams &french)
{
    const Real glass_h = frenchReducedGlassHeight(french);
    const Real hydro_speed = std::sqrt(std::max(cli.gravity_g * glass_h, TinyReal));
    if (cli.c0 <= TinyReal)
    {
        cli.c0 = Real(10) * hydro_speed;
    }
    if (cli.u_ref <= TinyReal)
    {
        cli.u_ref = hydro_speed;
    }
    if (cli.wall_rim_height < Real(0))
    {
        cli.wall_rim_height = Real(6) * french.dp;
    }
}

struct GlassEscapeCensus
{
    size_t escaped = 0;
    size_t total = 0;
    Real r_max = 0.0;
    Real z_min = 0.0;
    Real z_max = 0.0;
};

inline GlassEscapeCensus hostCountEscapedGlassParticles(BaseParticles &particles, const Vecd &center, Real radius,
                                                        Real z_bottom, Real z_top, Real dp)
{
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    GlassEscapeCensus census;
    census.total = particles.TotalRealParticles();
    if (census.total == 0)
    {
        return census;
    }
    census.z_min = pos[0][2];
    census.z_max = pos[0][2];
    const Real r_limit = radius + Real(2) * dp;
    const Real z_lo = z_bottom - Real(2) * dp;
    const Real z_hi = z_top + Real(2) * dp;
    for (size_t i = 0; i < census.total; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        census.r_max = std::max(census.r_max, r);
        census.z_min = std::min(census.z_min, pos[i][2]);
        census.z_max = std::max(census.z_max, pos[i][2]);
        if (r > r_limit || pos[i][2] < z_lo || pos[i][2] > z_hi)
        {
            ++census.escaped;
        }
    }
    return census;
}

inline void hostApplyHydrostaticDensity(BaseParticles &particles, Real rho0, Real c0, Real gravity_g, Real z_top)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "Density");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real inv_c2 = Real(1) / (c0 * c0 + TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        const Real depth = std::max(Real(0), z_top - pos[i][2]);
        rho[i] = rho0 * (Real(1) + gravity_g * depth * inv_c2);
    }
    syncVariableToDevice<Real>(particles, "Density");
}

inline void hostZeroVelocity(BaseParticles &particles)
{
    syncVariableToHost<Vecd>(particles, "Velocity");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        vel[i] = Vecd::Zero();
    }
    syncVariableToDevice<Vecd>(particles, "Velocity");
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

inline void hostSampleEulerQOntoParticles(BaseParticles &particles, const rh200::Rh200ScalarEulerianGrid &grid)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "JouleHeat");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *q = particles.getVariableDataByName<Real>("JouleHeat");
    const size_t n = particles.TotalRealParticles();
    for (size_t i = 0; i < n; ++i)
    {
        q[i] = grid.sampleHost(pos[i]);
    }
    syncVariableToDevice<Real>(particles, "JouleHeat");
}

inline Real hostThermalEnergy(BaseParticles &particles, Real rho, Real cp)
{
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *t = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    Real e = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        e += rho * cp * t[i] * vol[i];
    }
    return e;
}

inline bool depositEmQAndRecord(BaseParticles &em_particles, const OphelieGlassFieldNames &names,
                                const OphelieParameters &params, const FrenchCylindricalFrame &frame,
                                const FrenchLiteratureParameters &lit, rh200::Rh200ScalarEulerianGrid &q_grid,
                                ThermoEMCouplingState &state, Real p_em_recon)
{
    const std::string q_field = ophelieJouleHeatSourceFieldForThermal(names, params);
    syncVariableToHost<Real>(em_particles, q_field);
    syncVariableToHost<Vecd>(em_particles, "Position");
    syncVariableToHost<Real>(em_particles, "VolumetricMeasure");
    const Real *q = em_particles.getVariableDataByName<Real>(q_field);
    const Vecd *pos = em_particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = em_particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = em_particles.TotalRealParticles();

    state.q_rz.frame = frame;
    azimuthalDeposit(state.q_rz, pos, q, vol, n);
    fillEmptyAxisymmetricBins(state.q_rz, Real(0));

    StdVec<Vecd> em_pos(pos, pos + n);
    StdVec<Real> em_q(q, q + n);
    StdVec<Real> em_vol(vol, vol + n);
    q_grid.resetAccumulators();
    q_grid.depositScalarCloudInCell(em_pos, em_q, em_vol);
    q_grid.finalizeFromAccumulators();

    FrenchPowerMapBook book;
    book.p_em_recon = p_em_recon;
    book.p_axisymmetric_grid = integrateAxisymmetricPower(state.q_rz);
    book.p_euler_sample = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
    book.p_sph_particles = hostParticlePower(q, vol, n);
    book.p_euler_before_remap = book.p_euler_sample;
    if (lit.q_conservative_remap && book.p_euler_sample > TinyReal && p_em_recon > TinyReal)
    {
        book.remap_scale = p_em_recon / book.p_euler_sample;
        q_grid.scaleField(book.remap_scale);
        book.conservative_remap_applied = true;
        book.p_euler_sample = rh200::hostSamplePowerFromScalarGrid(q_grid, em_pos, em_vol);
        book.p_euler_after_remap = book.p_euler_sample;
        std::cout << "[ophelie][q-map] conservative remapping (not EM calibration): P_before="
                  << book.p_euler_before_remap << " scale=" << book.remap_scale << " P_after="
                  << book.p_euler_after_remap << std::endl;
    }
    logFrenchPowerMapBook(book);
    state.last_power = book;
    state.reconstructed_glass_power_w = p_em_recon;
    return true;
}

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
    if (local_cli.literature.coupling == ThermoEMCouplingMode::Periodic)
    {
        local_cli.enable_self_induction = local_cli.literature.enable_aind;
        local_cli.literature.target_glass_absorbed_power_w = local_cli.target_power_w;
    }
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
    OphelieGlassFieldNames glass_em_names;
    ThermoEMCouplingState coupling_state;
    rh200::Rh200ScalarEulerianGrid q_grid;
    const bool periodic_coupling = local_cli.literature.coupling == ThermoEMCouplingMode::Periodic;
    fs::create_directories("output");
    logFrenchLiteratureParameters(local_cli.literature);
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
        glass_em_names = glass_names;
        register_glass_em = makeUnique<RegisterOphelieGlassFields>(*glass_em, glass_em_names);
        glass_em_inner = makeUnique<Inner<>>(*glass_em);
        StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(*glass_em, glass_em_names,
                                                                                    params.sigma_glass_);
        assign_sigma.exec();

        OphelieFrenchEmJouleHeatOneWayResult probe;
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
        if (local_cli.enable_target_power)
        {
            (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
            runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
        }
        p_joule = probe.joule_power_w;
        phi_eq = probe.phi_eq_res_vol;
        coupling_state.calibrated_current_per_loop = french.coil.current_per_loop;
        coupling_state.coil_current_per_loop = french.coil.current_per_loop;
        coupling_state.current_frozen = local_cli.literature.em_control == EMControlMode::FixedCoilCurrent;
        coupling_state.generator_power_w = local_cli.literature.generator_power_w;
        coupling_state.target_glass_absorbed_power_w = local_cli.target_power_w;
        coupling_state.reconstructed_glass_power_w = p_joule;
        coupling_state.last_phi_eq_res = phi_eq;

        BaseParticles &particles = glass_em->getBaseParticles();
        const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_em_names, params);
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
                  << " generator_power=" << coupling_state.generator_power_w
                  << " target_glass_absorbed_power=" << coupling_state.target_glass_absorbed_power_w
                  << " reconstructed_glass_power=" << p_joule
                  << " coil_current=" << coupling_state.coil_current_per_loop << std::endl;

        if (periodic_coupling)
        {
            const Vecd &gc_em = french.glass_center;
            const FrenchCylindricalFrame frame = makeFrenchCylinderFrameFromGlass(
                gc_em, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
            const BoundingBoxd melt_bounds(
                Vecd(gc_em[0] - french.glass_radius, gc_em[1] - french.glass_radius,
                     gc_em[2] - french.glass_half_height),
                Vecd(gc_em[0] + french.glass_radius, gc_em[1] + french.glass_radius,
                     gc_em[2] + french.glass_half_height));
            q_grid.spec_ = rh200::makeRh200JouleHeatGridSpecFromBounds(melt_bounds, french.dp, Real(2));
            coupling_state.sigma_rz.frame = frame;
            coupling_state.q_rz.frame = frame;
            depositEmQAndRecord(particles, glass_em_names, params, frame, local_cli.literature, q_grid, coupling_state,
                                p_joule);
            syncVariableToHost<Real>(particles, glass_em_names.sigma);
            const Real *sigma0 = particles.getVariableDataByName<Real>(glass_em_names.sigma);
            const Vecd *pos0 = particles.getVariableDataByName<Vecd>("Position");
            const Real *vol0 = particles.getVariableDataByName<Real>("VolumetricMeasure");
            azimuthalDeposit(coupling_state.sigma_rz, pos0, sigma0, vol0, n);
            fillEmptyAxisymmetricBins(coupling_state.sigma_rz, params.sigma_glass_);
            coupling_state.sigma_em.assign(sigma0, sigma0 + n);
            sigmaFieldStats(coupling_state.sigma_em, vol0, n, coupling_state.last_sigma_min,
                            coupling_state.last_sigma_max, coupling_state.last_sigma_mean);
            coupling_state.last_update_time = 0.0;
            coupling_state.update_count = 1;
            writeAxisymmetricCsv("output/french_sigma_rz.csv", coupling_state.sigma_rz, "sigma_Sm");
            writeAxisymmetricCsv("output/french_q_rz.csv", coupling_state.q_rz, "Q_Wpm3");
            writeThermoEMCouplingCsvHeader("output/french_thermo_em_coupling.csv");
            appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, Real(0),
                                      local_cli.literature.em_control);
            writeFrenchEnergyBudgetCsvHeader("output/french_energy_budget.csv");
            std::cout << "[ophelie][thermo-em] initial coupling book written; mode="
                      << emControlModeName(local_cli.literature.em_control) << std::endl;
        }
    }

    //----------------------------------------------------------------------
    // Phase B — FluidBody WCSPH + Boussinesq + Natural thermal BC.
    // Domain is glass±wall only (not coil): Lattice scans the full SPHSystem
    // AABB, so a coil-sized box makes wall generation extremely slow.
    //----------------------------------------------------------------------
    finalizeNaturalHydroDefaults(local_cli, french);
    const Real wall_thickness = local_cli.wall_thickness_factor * french.dp;
    const Real rim_height = local_cli.wall_rim_height;
    const Real flow_margin = Real(4.0) * french.dp + wall_thickness;
    const Vecd &gc = french.glass_center;
    const Real flow_r = french.glass_radius + wall_thickness + flow_margin;
    const Real flow_down = french.glass_half_height + wall_thickness + flow_margin;
    const Real flow_up = french.glass_half_height + rim_height + wall_thickness + flow_margin;
    const BoundingBoxd flow_bounds(Vecd(gc[0] - flow_r, gc[1] - flow_r, gc[2] - flow_down),
                                   Vecd(gc[0] + flow_r, gc[1] + flow_r, gc[2] + flow_up));
    SPHSystem sph_system(flow_bounds, french.dp);
    sph_system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

    const OphelieFrenchGlassMaterialLaws laws = local_cli.literature.material;
    const Real t0 = local_cli.literature.t_initial_k;
    const Real rho0 = laws.rho.evaluate(t0);
    const Real mu = laws.mu.evaluate(t0);
    const Real cp = laws.cp.evaluate(t0);
    const Real k_th = laws.k.evaluate(t0);
    const Real beta = local_cli.literature.enable_boussinesq
                          ? (local_cli.beta > TinyReal ? local_cli.beta : laws.beta.evaluate(t0))
                          : Real(0);

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
                                   "WallBoundary", french, wall_thickness, wall_mesh_resolution, rim_height));
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
        if (periodic_coupling)
        {
            hostSampleEulerQOntoParticles(glass_particles, q_grid);
        }
    }

    OphelieThermalDiffusionOneWayOptions thermal_bc;
    thermal_bc.enable_french_natural_bc = true;
    thermal_bc.enable_cold_wall_dirichlet = false;
    thermal_bc.enable_diffusion = local_cli.literature.enable_thermal_diffusion;
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
    auto &t_min_reduce =
        main_methods.addReduceDynamics<OphelieThermalMinTemperatureReduceCK>(glass, kOphelieTemperatureField);

    StateDynamics<MainExecutionPolicy, ApplyFrozenQToTemperatureCK> apply_q(glass, rho0, cp, t0);
    StateDynamics<MainExecutionPolicy, ApplyOphelieThermalFrenchNaturalBcCK> apply_natural_bc(
        glass, kOphelieThermalDeltaTField, kOphelieTemperatureField, t0, rho0, cp, thermal_bc.boundary_shell_thickness,
        thermal_bc);
    UniquePtr<InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>>> laplace_temperature;
    UniquePtr<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>> apply_diffusion;
    if (thermal_bc.enable_diffusion)
    {
        laplace_temperature = makeUnique<InteractionDynamicsCK<MainExecutionPolicy, OpheliePairwiseLaplaceCK<Inner<>>>>(
            glass_inner, kOphelieTemperatureField, kOphelieThermalConductivityField, kOphelieThermalLaplaceTField,
            thermal_bc.pair_weight_regularization);
        apply_diffusion = makeUnique<StateDynamics<MainExecutionPolicy, ApplyLaplaceDiffusionDtCK>>(glass, rho0, cp,
                                                                                                   k_th, t0);
    }

    auto &glass_state_recorder = main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(glass);
    glass_state_recorder.addToWrite<Vecd>(glass, "Velocity");
    glass_state_recorder.addToWrite<Real>(glass, "Pressure");
    glass_state_recorder.addToWrite<Real>(glass, kOphelieTemperatureField);
    glass_state_recorder.addToWrite<Real>(glass, "JouleHeat");
    const bool write_vtp = local_cli.state_record_interval > TinyReal;

    const FrenchCylindricalFrame cyl_frame = makeFrenchCylinderFrameFromGlass(
        french.glass_center, french.glass_radius, french.glass_half_height, local_cli.literature.vertical_axis);
    Real prev_energy = hostThermalEnergy(glass_particles, rho0, cp);
    Real prev_energy_time = 0.0;

    update_wall_cell_linked_list.exec();
    update_glass_cell_linked_list.exec();
    update_glass_relations.exec();
    density_regularization.exec();
    const Real z_bottom = gc[2] - french.glass_half_height;
    const Real z_free = gc[2] + french.glass_half_height;
    const Real z_rim = z_free + rim_height;
    hostApplyHydrostaticDensity(glass_particles, rho0, local_cli.c0, local_cli.gravity_g, z_free);
    fluid_advection_step_setup.exec();
    constant_gravity.exec();

    TimeStepper &time_stepper = sph_solver.getTimeStepper();
    auto &advection_step = time_stepper.addTriggerByInterval(advection_time_step.exec());
    auto &state_recording =
        time_stepper.addTriggerByInterval(std::max(local_cli.state_record_interval, Real(1.0e-6)));
    size_t advection_steps = 0;

    const Real glass_h = frenchReducedGlassHeight(french);
    const Real alpha_th = k_th / (rho0 * cp + TinyReal);
    std::cout << "[ophelie][stage4.1] hydro: c0=" << local_cli.c0 << " u_ref=" << local_cli.u_ref
              << " rim=" << rim_height << " wall_thickness=" << wall_thickness
              << " settle_s=" << local_cli.hydro_settle_time
              << " gH/c0^2=" << (local_cli.gravity_g * glass_h) / (local_cli.c0 * local_cli.c0 + TinyReal)
              << " k=" << k_th << " alpha=" << alpha_th
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0) << std::endl;

    std::ofstream monitor("output/french_natural_convection_monitor.csv");
    monitor << "advection_step,physical_time_s,U_max,T_min,T_mean,T_max,dT,U_buoy,escaped,escaped_frac,"
               "P_joule,dE_dt,energy_residual,energy_residual_rel,total_heat_loss,em_updates,"
               "U_rms,U_theta_rms,U_z_rms,T_std,T_outer,T_center,T_bottom,T_top\n";
    writeFrenchSpatialStatsCsvHeader("output/french_spatial_stats.csv");

    while (time_stepper.getPhysicalTime() + TinyReal < local_cli.hydro_settle_time)
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);
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
        }
    }
    if (local_cli.hydro_settle_time > TinyReal)
    {
        const Real u_settle = u_max_reduce.exec();
        const GlassEscapeCensus settle_escape = hostCountEscapedGlassParticles(
            glass_particles, gc, french.glass_radius, z_bottom, z_rim, french.dp);
        std::cout << "[ophelie][stage4.1] hydro settle done: t=" << time_stepper.getPhysicalTime()
                  << " U_max=" << u_settle << " escaped=" << settle_escape.escaped << "/" << settle_escape.total
                  << " r_max=" << settle_escape.r_max << " z=[" << settle_escape.z_min << "," << settle_escape.z_max
                  << "]" << std::endl;
        hostZeroVelocity(glass_particles);
        time_stepper.setPhysicalTime(0.0);
        advection_steps = 0;
    }
    if (local_cli.literature.enable_boussinesq)
    {
        boussinesq_force.exec();
    }
    if (write_vtp)
    {
        glass_state_recorder.writeToFile(0);
    }

    std::cout << "[ophelie][stage4.1] flow start: end_time=" << local_cli.end_time
              << " n_glass=" << glass_particles.TotalRealParticles()
              << " n_wall=" << wall.getBaseParticles().TotalRealParticles() << " beta=" << beta << " mu=" << mu
              << " k=" << k_th << " c0=" << local_cli.c0 << " T0=" << t0
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0) << std::endl;

    while (!time_stepper.isEndTime(local_cli.end_time))
    {
        const Real acoustic_dt = time_stepper.incrementPhysicalTime(acoustic_time_step);
        acoustic_step_1st_half.exec(acoustic_dt);
        acoustic_step_2nd_half.exec(acoustic_dt);

        if (!local_cli.literature.zero_q)
        {
            apply_q.exec(acoustic_dt);
        }
        if (laplace_temperature && apply_diffusion)
        {
            laplace_temperature->exec();
            apply_diffusion->exec(acoustic_dt);
        }
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
            if (periodic_coupling)
            {
                hostSampleEulerQOntoParticles(glass_particles, q_grid);
                const Real t_now = time_stepper.getPhysicalTime();
                if (shouldUpdateEM(t_now, local_cli.literature.em_update_interval_s, coupling_state.last_update_time))
                {
                    syncVariableToHost<Real>(glass_particles, kOphelieTemperatureField);
                    syncVariableToHost<Vecd>(glass_particles, "Position");
                    syncVariableToHost<Real>(glass_particles, "VolumetricMeasure");
                    const Real *t_host = glass_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
                    const Vecd *pos_f = glass_particles.getVariableDataByName<Vecd>("Position");
                    const Real *vol_f = glass_particles.getVariableDataByName<Real>("VolumetricMeasure");
                    const size_t n_f = glass_particles.TotalRealParticles();
                    BaseParticles &em_particles = glass_em->getBaseParticles();
                    syncVariableToHost<Vecd>(em_particles, "Position");
                    syncVariableToHost<Real>(em_particles, "VolumetricMeasure");
                    const Vecd *pos_em = em_particles.getVariableDataByName<Vecd>("Position");
                    const Real *vol_em = em_particles.getVariableDataByName<Real>("VolumetricMeasure");
                    const size_t n_em = em_particles.TotalRealParticles();
                    StdVec<Real> sigma_on_em;
                    std::string sigma_error;
                    if (!mapFluidTemperatureToUnderRelaxedEmSigma(local_cli.literature, cyl_frame, pos_f, vol_f,
                                                                  t_host, n_f, pos_em, vol_em, n_em, coupling_state,
                                                                  sigma_on_em, sigma_error))
                    {
                        std::cerr << "[ophelie][thermo-em] abort: " << sigma_error << std::endl;
                        return 1;
                    }
                    hostAssignSigmaField(em_particles, glass_em_names.sigma, sigma_on_em);

                    OphelieFrenchEmJouleHeatOneWayResult probe;
                    runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                        *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
                    if (local_cli.literature.em_control == EMControlMode::FixedAbsorbedPower)
                    {
                        (void)calibrateFrenchCoilCurrentToTargetPower(french, params, probe.joule_power_w);
                        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
                            *glass_em, *glass_em_inner, glass_em_names, params, french, probe);
                    }
                    p_joule = probe.joule_power_w;
                    phi_eq = probe.phi_eq_res_vol;
                    coupling_state.coil_current_per_loop = french.coil.current_per_loop;
                    coupling_state.last_phi_eq_res = phi_eq;
                    depositEmQAndRecord(em_particles, glass_em_names, params, cyl_frame, local_cli.literature, q_grid,
                                        coupling_state, p_joule);
                    hostSampleEulerQOntoParticles(glass_particles, q_grid);
                    coupling_state.last_update_time = t_now;
                    ++coupling_state.update_count;
                    writeAxisymmetricCsv("output/french_sigma_rz.csv", coupling_state.sigma_rz, "sigma_Sm");
                    writeAxisymmetricCsv("output/french_q_rz.csv", coupling_state.q_rz, "Q_Wpm3");
                    appendThermoEMCouplingCsv("output/french_thermo_em_coupling.csv", coupling_state, t_now,
                                              local_cli.literature.em_control);
                    std::cout << "[ophelie][thermo-em] update=" << coupling_state.update_count << " t=" << t_now
                              << " reconstructed_glass_power=" << p_joule
                              << " coil_current=" << coupling_state.coil_current_per_loop
                              << " phi_eq_res_vol=" << phi_eq << std::endl;
                }
            }
            if (write_vtp && state_recording())
            {
                glass_state_recorder.writeToFile(advection_steps);
            }
            if (advection_steps % 20 == 0)
            {
                const Real t_now = time_stepper.getPhysicalTime();
                const Real energy = hostThermalEnergy(glass_particles, rho0, cp);
                const Real dt_e = t_now - prev_energy_time;
                const Real dE_dt = dt_e > TinyReal ? (energy - prev_energy) / dt_e : Real(0);
                Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
                hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc,
                                                              thermal_bc.boundary_shell_thickness, side_w, bottom_w,
                                                              free_conv_w, free_rad_w, total_loss_w);
                const Real u_now = u_max_reduce.exec();
                const Real t_mean_now = t_mean_reduce.exec();
                const Real t_max_now = t_max_reduce.exec();
                const Real t_min_now = t_min_reduce.exec();
                const Real dT_now = t_max_now - t_min_now;
                const Real u_buoy = frenchBuoyancySpeedScale(local_cli.gravity_g, beta, dT_now, glass_h);
                const Real residual =
                    frenchEnergyResidual(dE_dt, p_joule, side_w, free_conv_w + free_rad_w, bottom_w);
                const Real residual_rel = residual / (std::abs(p_joule) + TinyReal);
                const FrenchMeltSpatialStats spatial =
                    hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
                appendFrenchSpatialStatsCsv("output/french_spatial_stats.csv", t_now, spatial);
                writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame,
                                           kOphelieTemperatureField);
                if (periodic_coupling)
                {
                    appendFrenchEnergyBudgetCsv("output/french_energy_budget.csv", t_now, energy, dE_dt, p_joule,
                                                side_w, free_conv_w + free_rad_w, bottom_w, t_min_now, t_mean_now,
                                                t_max_now, u_now, u_buoy,
                                                thermal_bc.enable_diffusion ? 1 : 0);
                }
                prev_energy = energy;
                prev_energy_time = t_now;
                const GlassEscapeCensus mid_escape = hostCountEscapedGlassParticles(
                    glass_particles, gc, french.glass_radius, z_bottom, z_rim, french.dp);
                const Real mid_frac =
                    mid_escape.total > 0 ? static_cast<Real>(mid_escape.escaped) / static_cast<Real>(mid_escape.total)
                                         : Real(0);
                monitor << advection_steps << "," << t_now << "," << u_now << "," << t_min_now << "," << t_mean_now
                        << "," << t_max_now << "," << dT_now << "," << u_buoy << "," << mid_escape.escaped << ","
                        << mid_frac << "," << p_joule << "," << dE_dt << "," << residual << "," << residual_rel << ","
                        << total_loss_w << "," << coupling_state.update_count << "," << spatial.u_rms << ","
                        << spatial.u_theta_rms << "," << spatial.u_z_rms << "," << spatial.t_std << ","
                        << spatial.t_outer << "," << spatial.t_center << "," << spatial.t_bottom << ","
                        << spatial.t_top << "\n";
                monitor.flush();
                std::cout << "[ophelie][stage4.1] N=" << advection_steps << " t=" << t_now << " U_max=" << u_now
                          << " U_rms=" << spatial.u_rms << " U_th=" << spatial.u_theta_rms << " U_z=" << spatial.u_z_rms
                          << " U_buoy=" << u_buoy << " T_min=" << t_min_now << " T_mean=" << t_mean_now
                          << " T_max=" << t_max_now << " dT=" << dT_now << " T_std=" << spatial.t_std
                          << " T_out=" << spatial.t_outer << " T_ctr=" << spatial.t_center
                          << " escaped=" << mid_escape.escaped << "/" << mid_escape.total
                          << " residual=" << residual << " residual_rel=" << residual_rel
                          << " em_updates=" << coupling_state.update_count << std::endl;
            }
        }
    }

    const Real u_max = u_max_reduce.exec();
    const Real t_mean = t_mean_reduce.exec();
    const Real t_max = t_max_reduce.exec();
    const Real t_min = t_min_reduce.exec();
    const Real dT = t_max - t_min;
    const Real u_buoy = frenchBuoyancySpeedScale(local_cli.gravity_g, beta, dT, glass_h);
    const FrenchMeltSpatialStats spatial_end =
        hostFrenchMeltSpatialStats(glass_particles, cyl_frame, kOphelieTemperatureField);
    writeFluidTemperatureRzCsv("output/french_T_rz.csv", glass_particles, cyl_frame, kOphelieTemperatureField);
    const GlassEscapeCensus escape = hostCountEscapedGlassParticles(glass_particles, gc, french.glass_radius, z_bottom,
                                                                   z_rim, french.dp);
    const Real escaped_frac =
        escape.total > 0 ? static_cast<Real>(escape.escaped) / static_cast<Real>(escape.total) : Real(0);
    Real side_w = 0, bottom_w = 0, free_conv_w = 0, free_rad_w = 0, total_loss_w = 0;
    hostOphelieThermalFrenchNaturalHeatLossPowers(glass_particles, thermal_bc, thermal_bc.boundary_shell_thickness,
                                                  side_w, bottom_w, free_conv_w, free_rad_w, total_loss_w);
    if (write_vtp)
    {
        glass_state_recorder.writeToFile(advection_steps + 1);
    }
    monitor.close();

    const Real power_rel_err =
        std::abs(p_joule - local_cli.target_power_w) / (local_cli.target_power_w + TinyReal);
    const bool em_ok =
        std::isfinite(p_joule) &&
        (periodic_coupling ? std::isfinite(phi_eq) : (power_rel_err < Real(1.0e-2)));
    const bool flow_ok = std::isfinite(u_max) && std::isfinite(t_mean) && std::isfinite(t_max) &&
                         u_max <= local_cli.c0 && escaped_frac <= local_cli.escaped_frac_max;
    const bool thermal_ok = total_loss_w > TinyReal && side_w > TinyReal && bottom_w > TinyReal &&
                            (free_conv_w + free_rad_w) > TinyReal;
    const bool buoyancy_ok = true;
    const bool coupling_ok =
        !periodic_coupling || (coupling_state.update_count >= 1 && std::isfinite(coupling_state.last_phi_eq_res));
    const bool passed = em_ok && flow_ok && thermal_ok && buoyancy_ok && coupling_ok;

    std::cout << "test_3d_ophelie_french_natural_convection_frozen_q"
              << " P_joule_W=" << p_joule << " power_rel_err=" << power_rel_err << " phi_eq_res_vol=" << phi_eq
              << " end_time=" << local_cli.end_time << " U_max=" << u_max << " U_rms=" << spatial_end.u_rms
              << " U_th=" << spatial_end.u_theta_rms << " U_z=" << spatial_end.u_z_rms << " U_buoy=" << u_buoy
              << " T_min=" << t_min << " T_mean=" << t_mean << " T_max=" << t_max << " dT=" << dT
              << " T_std=" << spatial_end.t_std << " T_out=" << spatial_end.t_outer << " T_ctr=" << spatial_end.t_center
              << " escaped=" << escape.escaped << "/" << escape.total << " escaped_frac=" << escaped_frac
              << " r_max=" << escape.r_max << " z=[" << escape.z_min << "," << escape.z_max << "]"
              << " wall_loss_side_W=" << side_w << " wall_loss_bottom_W=" << bottom_w
              << " free_conv_loss_W=" << free_conv_w << " free_rad_loss_W=" << free_rad_w
              << " total_heat_loss_W=" << total_loss_w
              << " thermal_diffusion=" << (thermal_bc.enable_diffusion ? 1 : 0) << " em_ok=" << (em_ok ? 1 : 0)
              << " flow_ok=" << (flow_ok ? 1 : 0) << " thermal_ok=" << (thermal_ok ? 1 : 0)
              << " buoyancy_ok=" << (buoyancy_ok ? 1 : 0) << " coupling_updates=" << coupling_state.update_count
              << " coupling=" << thermoEMCouplingModeName(local_cli.literature.coupling)
              << " em_control=" << emControlModeName(local_cli.literature.em_control)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
