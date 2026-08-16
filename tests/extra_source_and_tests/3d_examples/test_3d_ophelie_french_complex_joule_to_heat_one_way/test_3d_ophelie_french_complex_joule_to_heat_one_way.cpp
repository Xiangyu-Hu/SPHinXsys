/**
 * @file test_3d_ophelie_french_complex_joule_to_heat_one_way.cpp
 * @brief French reload: complex edge-flux EM (optional A_glass Picard) → JouleHeat → one-way thermal.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_sigma_t_coupling.h"
#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"
#include "electromagnetic_ophelie_thermal_vtp.h"
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

struct ThermalLocalCli
{
    bool use_reload = false;
    std::string reload_dir;
    Real thermal_dt = 1.0;
    size_t thermal_steps = 3;
    bool thermal_diffusion = false;
    /** dirichlet_regression (default with diffusion) | french_natural_production */
    std::string thermal_bc_name = "dirichlet_regression";
    Real h_side = 300.0;
    Real h_bottom = 35.0;
    Real h_free = 20.0;
    Real emissivity = 0.8;
    Real t_cool = 300.0;
    Real t_ambient = 300.0;
    Real t_rad_ambient = 300.0;
    OphelieThermalMaterialPreset material_preset = OphelieThermalMaterialPreset::ReducedPrototype;
    bool t0_user_set = false;
    bool rho_user_set = false;
    bool cp_user_set = false;
    bool k_user_set = false;
    Real t0_override = 0.0;
    Real rho_override = 0.0;
    Real cp_override = 0.0;
    Real k_override = 0.0;
    bool thermal_state_recording = false;
    size_t thermal_record_interval = 0;
    bool enable_sigma_t = false;
    size_t sigma_t_max_iter = 8;
    Real sigma_t_relax = 0.3;
    Real sigma_t_tol = 0.05;
    Real sigma_t_seed_delta_k = 100.0;
    /** provisional | french-natural | thesis-iii12 | thesis-iv10 */
    std::string sigma_law_name = "provisional";
    bool enable_target_power = false;
    Real target_power_w = 50000.0;
};

inline OphelieThermalMaterialPreset parseThermalMaterialPreset(const char *value)
{
    if (std::strcmp(value, "literature") == 0 || std::strcmp(value, "jacoutot") == 0 ||
        std::strcmp(value, "jacoutot_table1_1473K") == 0)
    {
        return OphelieThermalMaterialPreset::JacoutotTable1_1473K;
    }
    return OphelieThermalMaterialPreset::ReducedPrototype;
}

inline void applyThermalMaterialFromCli(const ThermalLocalCli &cli, OphelieJouleHeatOneWayMaterialProps &material)
{
    applyOphelieThermalMaterialPreset(material, cli.material_preset);
    if (cli.t0_user_set)
    {
        material.t_initial = cli.t0_override;
    }
    if (cli.rho_user_set)
    {
        material.rho = cli.rho_override;
    }
    if (cli.cp_user_set)
    {
        material.cp = cli.cp_override;
    }
    if (cli.k_user_set)
    {
        material.k = cli.k_override;
    }
}

inline OphelieTemperatureLaw selectSigmaTemperatureLaw(const std::string &name, bool &paper_digitized)
{
    if (name == "french-natural" || name == "cep2008" || name == "cep2008-nat" || name == "natural")
    {
        paper_digitized = true;
        return makeFrenchNaturalJournalAnchoredSigmaTemperatureLaw();
    }
    if (name == "thesis-iii12" || name == "iii12" || name == "fig-iii12")
    {
        paper_digitized = true;
        return makeJacoutotThesisSigmaTemperatureLawIII12();
    }
    if (name == "thesis-iv10" || name == "iv10" || name == "fig-iv10")
    {
        paper_digitized = true;
        return makeJacoutotThesisSigmaTemperatureLawIV10();
    }
    paper_digitized = false;
    return makeProvisionalJacoutotLikeSigmaTemperatureLaw();
}

inline void applyThermalLocalCli(int ac, char *av[], OphelieFrenchReducedCaseParams &french, ThermalLocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            cli.use_reload = true;
            french.dp = 0.02;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            cli.reload_dir = std::string(av[i] + 13);
            cli.use_reload = true;
        }
        else if (std::strncmp(av[i], "--thermal-dt=", 13) == 0)
        {
            cli.thermal_dt = static_cast<Real>(std::atof(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--thermal-steps=", 16) == 0)
        {
            cli.thermal_steps = static_cast<size_t>(std::atoi(av[i] + 16));
        }
        else if (std::strcmp(av[i], "--thermal-diffusion=1") == 0 || std::strcmp(av[i], "--thermal-diffusion") == 0)
        {
            cli.thermal_diffusion = true;
        }
        else if (std::strncmp(av[i], "--thermal-bc=", 13) == 0)
        {
            cli.thermal_bc_name = std::string(av[i] + 13);
            cli.thermal_diffusion = true;
        }
        else if (std::strncmp(av[i], "--h-side=", 9) == 0)
        {
            cli.h_side = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--h-bottom=", 11) == 0)
        {
            cli.h_bottom = static_cast<Real>(std::atof(av[i] + 11));
        }
        else if (std::strncmp(av[i], "--h-free=", 9) == 0)
        {
            cli.h_free = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--emissivity=", 13) == 0)
        {
            cli.emissivity = static_cast<Real>(std::atof(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--t-cool=", 9) == 0)
        {
            cli.t_cool = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--t-ambient=", 12) == 0)
        {
            cli.t_ambient = static_cast<Real>(std::atof(av[i] + 12));
        }
        else if (std::strncmp(av[i], "--t-rad-ambient=", 16) == 0)
        {
            cli.t_rad_ambient = static_cast<Real>(std::atof(av[i] + 16));
        }
        else if (std::strncmp(av[i], "--thermal-material=", 19) == 0)
        {
            cli.material_preset = parseThermalMaterialPreset(av[i] + 19);
        }
        else if (std::strcmp(av[i], "--use-literature-thermal=1") == 0 ||
                 std::strcmp(av[i], "--use-literature-thermal") == 0)
        {
            cli.material_preset = OphelieThermalMaterialPreset::JacoutotTable1_1473K;
        }
        else if (std::strncmp(av[i], "--thermal-t0=", 13) == 0)
        {
            cli.t0_user_set = true;
            cli.t0_override = static_cast<Real>(std::atof(av[i] + 13));
        }
        else if (std::strncmp(av[i], "--rho=", 6) == 0)
        {
            cli.rho_user_set = true;
            cli.rho_override = static_cast<Real>(std::atof(av[i] + 6));
        }
        else if (std::strncmp(av[i], "--cp=", 5) == 0)
        {
            cli.cp_user_set = true;
            cli.cp_override = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--k=", 4) == 0)
        {
            cli.k_user_set = true;
            cli.k_override = static_cast<Real>(std::atof(av[i] + 4));
        }
        else if (std::strcmp(av[i], "--thermal-state-recording=1") == 0 ||
                 std::strcmp(av[i], "--thermal-state-recording") == 0 ||
                 std::strcmp(av[i], "--state-recording=1") == 0)
        {
            cli.thermal_state_recording = true;
        }
        else if (std::strncmp(av[i], "--thermal-record-interval=", 26) == 0)
        {
            cli.thermal_record_interval = static_cast<size_t>(std::atoi(av[i] + 26));
        }
        else if (std::strcmp(av[i], "--sigma-t=1") == 0 || std::strcmp(av[i], "--sigma-t") == 0)
        {
            cli.enable_sigma_t = true;
        }
        else if (std::strncmp(av[i], "--sigma-t-max-iter=", 19) == 0)
        {
            cli.sigma_t_max_iter = static_cast<size_t>(std::atoi(av[i] + 19));
        }
        else if (std::strncmp(av[i], "--sigma-t-relax=", 16) == 0)
        {
            cli.sigma_t_relax = static_cast<Real>(std::atof(av[i] + 16));
        }
        else if (std::strncmp(av[i], "--sigma-t-tol=", 14) == 0)
        {
            cli.sigma_t_tol = static_cast<Real>(std::atof(av[i] + 14));
        }
        else if (std::strncmp(av[i], "--sigma-t-seed-delta-k=", 23) == 0)
        {
            cli.sigma_t_seed_delta_k = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strncmp(av[i], "--sigma-law=", 12) == 0)
        {
            cli.sigma_law_name = std::string(av[i] + 12);
        }
        else if (std::strncmp(av[i], "--target-power=", 15) == 0)
        {
            cli.enable_target_power = true;
            cli.target_power_w = static_cast<Real>(std::atof(av[i] + 15));
        }
        else if (std::strcmp(av[i], "--target-power") == 0)
        {
            cli.enable_target_power = true;
            cli.target_power_w = 50000.0;
        }
    }
    if (cli.use_reload && cli.reload_dir.empty())
    {
        cli.reload_dir = resolveDefaultFrenchReloadFolder();
    }
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
    ThermalLocalCli local_cli;
    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);
    applyThermalLocalCli(ac, av, french, local_cli);

    const StdVec<std::string> french_filtered = filterFrenchReducedCommandLine(ac, av, french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    OphelieTestCliOptions cli_options;
    StdVec<char *> french_av;
    french_av.reserve(french_filtered.size());
    for (auto &argument : french_filtered)
    {
        french_av.push_back(const_cast<char *>(argument.c_str()));
    }
    const int french_ac = static_cast<int>(french_av.size());
    (void)filterOphelieTestCommandLine(french_ac, french_av.data(), params, cli_options);

    if (!cli_options.reload_dir.empty())
    {
        local_cli.reload_dir = cli_options.reload_dir;
        local_cli.use_reload = true;
    }

    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);
    if (local_cli.enable_target_power)
    {
        // Stage 3.5: calibrate coil current to P_target; do not use field power scaling.
        params.target_joule_power_ = local_cli.target_power_w;
        params.enable_power_scaling_ = false;
    }
    logOphelieFinalParams(params, cli_options);

    if (!local_cli.use_reload || !reloadXmlExists(local_cli.reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_complex_joule_to_heat_one_way: Reload.xml required under \""
                  << local_cli.reload_dir << "\"\n";
        return 1;
    }

    const BoundingBoxd bounds = frenchReducedDomainBounds(french, 3.0 * french.dp);
    SPHSystem sph_system(bounds, french.dp);
    sph_system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);
    glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    OphelieJouleHeatOneWayMaterialProps material;
    applyThermalMaterialFromCli(local_cli, material);
    if (local_cli.enable_sigma_t && !local_cli.t0_user_set &&
        local_cli.material_preset == OphelieThermalMaterialPreset::ReducedPrototype)
    {
        // sigma(T) Arrhenius is anchored near melt temperatures; default to Jacoutot T0.
        applyOphelieThermalMaterialPreset(material, OphelieThermalMaterialPreset::JacoutotTable1_1473K);
        local_cli.material_preset = OphelieThermalMaterialPreset::JacoutotTable1_1473K;
    }

    OphelieThermalDiffusionOneWayOptions thermal_options;
    thermal_options.enable_diffusion = local_cli.thermal_diffusion;
    const bool natural_bc = local_cli.thermal_bc_name == "french-natural" ||
                            local_cli.thermal_bc_name == "french_natural" ||
                            local_cli.thermal_bc_name == "natural" ||
                            local_cli.thermal_bc_name == "production";
    if (natural_bc)
    {
        thermal_options.enable_diffusion = true;
        thermal_options.enable_french_natural_bc = true;
        thermal_options.enable_cold_wall_dirichlet = false;
        thermal_options.h_side = local_cli.h_side;
        thermal_options.h_bottom = local_cli.h_bottom;
        thermal_options.h_free = local_cli.h_free;
        thermal_options.emissivity = local_cli.emissivity;
        thermal_options.t_cool = local_cli.t_cool;
        thermal_options.t_ambient = local_cli.t_ambient;
        thermal_options.t_rad_ambient = local_cli.t_rad_ambient;
        local_cli.thermal_diffusion = true;
        if (!local_cli.t0_user_set &&
            local_cli.material_preset == OphelieThermalMaterialPreset::ReducedPrototype)
        {
            applyOphelieThermalMaterialPreset(material, OphelieThermalMaterialPreset::JacoutotTable1_1473K);
            local_cli.material_preset = OphelieThermalMaterialPreset::JacoutotTable1_1473K;
        }
    }
    else
    {
        // Default with --thermal-diffusion: Dirichlet regression (Stage 3.3).
        thermal_options.enable_cold_wall_dirichlet = local_cli.thermal_diffusion;
        thermal_options.enable_french_natural_bc = false;
    }
    thermal_options.boundary_width_factor = params.phi_boundary_distance_factor_;

    OphelieThermalVtpRecordingOptions thermal_recording;
    thermal_recording.enabled = local_cli.thermal_state_recording;
    thermal_recording.sph_system = &sph_system;
    thermal_recording.glass_body = &glass_body;
    thermal_recording.names = &glass_names;
    thermal_recording.params = &params;
    thermal_recording.record_interval = local_cli.thermal_record_interval;
    thermal_recording.include_diffusion_fields = local_cli.thermal_diffusion;
    const OphelieThermalVtpRecordingOptions *recording_ptr =
        local_cli.thermal_state_recording ? &thermal_recording : nullptr;

    OphelieFrenchEmJouleHeatOneWayResult pipeline;
    OphelieSigmaTCouplingResult sigma_t;
    Real boundary_compliance = 1.0;
    Real max_temperature = material.t_initial;
    Real joule_power_raw_before_calibrate = 0.0;
    Real coil_current_scale_to_target = 1.0;
    if (local_cli.enable_target_power)
    {
        // Probe EM at reference current, scale I ~ sqrt(P_target/P_raw), then re-solve below.
        OphelieFrenchEmJouleHeatOneWayResult probe;
        runFrenchReducedEmOrSelfInductionForThermalHandoff<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, french, probe);
        joule_power_raw_before_calibrate = probe.joule_power_w;
        coil_current_scale_to_target =
            calibrateFrenchCoilCurrentToTargetPower(french, params, joule_power_raw_before_calibrate);
    }
    if (local_cli.enable_sigma_t)
    {
        // σ(T) outer loop keeps lumped heating for auditable σ/P convergence.
        // Optional --thermal-diffusion then verifies conduction + cold wall on frozen final Q.
        OphelieSigmaTCouplingOptions sigma_options;
        sigma_options.max_iterations = local_cli.sigma_t_max_iter;
        sigma_options.sigma_relaxation = local_cli.sigma_t_relax;
        sigma_options.relative_tol = local_cli.sigma_t_tol;
        sigma_options.thermal_dt = local_cli.thermal_dt;
        sigma_options.thermal_steps_per_em =
            local_cli.thermal_diffusion ? size_t(1) : local_cli.thermal_steps;
        sigma_options.temperature_seed_delta_k = local_cli.sigma_t_seed_delta_k;
        bool paper_digitized_sigma = false;
        const OphelieTemperatureLaw sigma_law =
            selectSigmaTemperatureLaw(local_cli.sigma_law_name, paper_digitized_sigma);
        sigma_t = runFrenchReducedSigmaTEmThermalCoupling<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, french, material, sigma_law, sigma_options,
            paper_digitized_sigma);
        pipeline = sigma_t.last_handoff;
        max_temperature = material.t_initial + pipeline.thermal.max_delta_t;
        if (local_cli.thermal_diffusion)
        {
            BaseParticles &particles = glass_body.getBaseParticles();
            const std::string q_field = ophelieJouleHeatSourceFieldForThermal(glass_names, params);
            const OphelieThermalDiffusionOneWayStepResult diffusion =
                applyOphelieFrozenQDiffusionFromUniformT0<MainExecutionPolicy>(
                    glass_body, *glass_inner, particles, q_field, french, local_cli.thermal_dt,
                    local_cli.thermal_steps, material, thermal_options, recording_ptr);
            copyOphelieJouleHeatOneWayStepResult(pipeline.thermal, diffusion);
            boundary_compliance = diffusion.boundary_dirichlet_compliance;
            max_temperature = diffusion.max_temperature;
        }
        else if (local_cli.thermal_state_recording)
        {
            writeOphelieThermalBodyStatesVtp(sph_system, glass_body, glass_names, params, local_cli.thermal_steps,
                                             false);
        }
    }
    else if (local_cli.thermal_diffusion)
    {
        pipeline = runFrenchReducedEmThenJouleHeatDiffusionOneWay<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, french, local_cli.thermal_dt, local_cli.thermal_steps,
            material, thermal_options, recording_ptr);
        boundary_compliance = pipeline.thermal.boundary_dirichlet_compliance;
        max_temperature = pipeline.thermal.max_temperature;
    }
    else
    {
        pipeline = runFrenchReducedEmThenJouleHeatOneWay<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, french, local_cli.thermal_dt, local_cli.thermal_steps,
            material);
        if (local_cli.thermal_state_recording)
        {
            writeOphelieThermalBodyStatesVtp(sph_system, glass_body, glass_names, params, local_cli.thermal_steps,
                                             false);
        }
    }

    const size_t n = glass_body.getBaseParticles().TotalRealParticles();
    // Dual-layer φ gate (Stage 3.4 closeout): preferred from CLI/params; engineering hard ≥ 2.5e-4.
    const Real em_phi_preferred =
        pipeline.used_self_induction ? ophelieSelfInductionPicardPhiEqResTolerance(params) : Real(0.01);
    const Real em_phi_hard = std::max(em_phi_preferred, Real(2.5e-4));
    const bool j_ok =
        !pipeline.used_self_induction ||
        (pipeline.final_j_rel_change < params.self_induction_j_tolerance_);
    const bool phi_preferred_ok = pipeline.phi_eq_res_vol < em_phi_preferred;
    const bool phi_ok = pipeline.phi_eq_res_vol < em_phi_hard;
    const bool phi_gate_warn = phi_ok && !phi_preferred_ok;
    const Real power_target_rel_err =
        local_cli.enable_target_power
            ? std::abs(pipeline.joule_power_w - local_cli.target_power_w) / (local_cli.target_power_w + TinyReal)
            : Real(0);
    const bool power_ok =
        !local_cli.enable_target_power || power_target_rel_err < Real(1.0e-3);
    const bool em_ok = n > 0 && std::isfinite(pipeline.joule_power_w) && pipeline.joule_power_w > TinyReal && phi_ok &&
                       j_ok && power_ok &&
                       (!pipeline.used_self_induction || pipeline.a_ind_over_a_coil > TinyReal);
    const Real vol_weighted_rel_err =
        std::abs(pipeline.thermal.vol_weighted_delta_t - pipeline.thermal.vol_weighted_expected_delta_t) /
        (pipeline.thermal.vol_weighted_expected_delta_t + TinyReal);
    const Real expected_energy_j =
        pipeline.joule_power_w * local_cli.thermal_dt * static_cast<Real>(local_cli.thermal_steps);
    const Real energy_vs_power_rel_err =
        std::abs(pipeline.thermal.total_thermal_energy_j - expected_energy_j) / (expected_energy_j + TinyReal);
    const Real energy_cap_rel_err =
        (pipeline.thermal.total_thermal_energy_j - pipeline.thermal.total_joule_energy_j) /
        (pipeline.thermal.total_joule_energy_j + TinyReal);

    bool thermal_ok = false;
    if (thermal_options.enable_french_natural_bc)
    {
        thermal_ok = std::isfinite(max_temperature) && std::isfinite(pipeline.thermal.total_heat_loss_w) &&
                     pipeline.thermal.total_heat_loss_w > TinyReal &&
                     pipeline.thermal.wall_heat_loss_side_w > TinyReal &&
                     pipeline.thermal.wall_heat_loss_bottom_w > TinyReal &&
                     (pipeline.thermal.free_surface_conv_loss_w + pipeline.thermal.free_surface_rad_loss_w) >
                         TinyReal &&
                     energy_cap_rel_err <= Real(1.0e-6);
    }
    else if (local_cli.thermal_diffusion)
    {
        thermal_ok = pipeline.thermal.max_delta_t > TinyReal && std::isfinite(max_temperature) &&
                     max_temperature > material.t_initial + TinyReal && boundary_compliance > Real(0.90) &&
                     energy_cap_rel_err <= Real(1.0e-6);
    }
    else
    {
        thermal_ok = pipeline.thermal.max_delta_t > TinyReal && vol_weighted_rel_err < Real(0.05) &&
                     energy_vs_power_rel_err < Real(0.05) && pipeline.thermal.energy_balance_rel_err < Real(0.05) &&
                     pipeline.thermal.closure_mismatch_vol_fraction < Real(0.01);
    }

    const bool sigma_spatial_ok =
        !local_cli.enable_sigma_t ||
        (sigma_t.sigma_max > sigma_t.sigma_min * Real(1.05) && std::isfinite(sigma_t.sigma_mean));
    const bool sigma_coupling_ok =
        !local_cli.enable_sigma_t || (sigma_t.converged && sigma_t.coupling_iterations > 0);
    const bool passed = em_ok && thermal_ok && sigma_spatial_ok && sigma_coupling_ok;
    if (phi_gate_warn)
    {
        std::cerr << "[ophelie][phi-gate] engineering pass: phi_eq_res_vol=" << pipeline.phi_eq_res_vol
                  << " preferred=" << em_phi_preferred << " hard=" << em_phi_hard << std::endl;
    }

    std::cout << "test_3d_ophelie_french_complex_joule_to_heat_one_way particles=reload n=" << n
              << " dp=" << french.dp << " edge_flux_complex=1"
              << " self_induction=" << (pipeline.used_self_induction ? 1 : 0)
              << " self_induction_iters=" << pipeline.self_induction_iterations
              << " final_J_rel=" << pipeline.final_j_rel_change << " j_ok=" << (j_ok ? 1 : 0)
              << " picard_converged=" << (pipeline.picard_converged ? 1 : 0)
              << " A_ind_over_A_coil=" << pipeline.a_ind_over_a_coil
              << " B_ind_over_B_coil=" << pipeline.b_ind_over_b_coil
              << " sigma_t=" << (local_cli.enable_sigma_t ? 1 : 0)
              << " sigma_t_iters=" << sigma_t.coupling_iterations
              << " sigma_t_converged=" << (sigma_t.converged ? 1 : 0)
              << " sigma_mean=" << sigma_t.sigma_mean << " sigma_min=" << sigma_t.sigma_min
              << " sigma_max=" << sigma_t.sigma_max << " sigma_rel=" << sigma_t.final_sigma_rel
              << " power_rel=" << sigma_t.final_power_rel
              << " paper_digitized_sigma=" << (sigma_t.paper_digitized_sigma_law ? 1 : 0)
              << " sigma_law=" << local_cli.sigma_law_name
              << " thermal_diffusion=" << (local_cli.thermal_diffusion ? 1 : 0)
              << " thermal_bc=" << ophelieThermalBoundaryModeName(ophelieThermalBoundaryModeFromOptions(thermal_options))
              << " h_side=" << thermal_options.h_side << " h_bottom=" << thermal_options.h_bottom
              << " h_free=" << thermal_options.h_free << " emissivity=" << thermal_options.emissivity
              << " t_cool=" << thermal_options.t_cool
              << " wall_loss_side_W=" << pipeline.thermal.wall_heat_loss_side_w
              << " wall_loss_bottom_W=" << pipeline.thermal.wall_heat_loss_bottom_w
              << " free_conv_loss_W=" << pipeline.thermal.free_surface_conv_loss_w
              << " free_rad_loss_W=" << pipeline.thermal.free_surface_rad_loss_w
              << " total_heat_loss_W=" << pipeline.thermal.total_heat_loss_w
              << " thermal_material=" << ophelieThermalMaterialPresetName(local_cli.material_preset)
              << " T0=" << material.t_initial << " rho=" << material.rho << " cp=" << material.cp << " k=" << material.k
              << " target_power=" << (local_cli.enable_target_power ? local_cli.target_power_w : Real(0))
              << " P_raw_before_calibrate=" << joule_power_raw_before_calibrate
              << " I_scale_to_target=" << coil_current_scale_to_target
              << " power_target_rel_err=" << power_target_rel_err << " power_ok=" << (power_ok ? 1 : 0)
              << " P_joule_W=" << pipeline.joule_power_w << " phi_eq_res_vol=" << pipeline.phi_eq_res_vol
              << " em_phi_preferred=" << em_phi_preferred << " em_phi_hard=" << em_phi_hard
              << " phi_preferred_ok=" << (phi_preferred_ok ? 1 : 0) << " phi_gate_warn=" << (phi_gate_warn ? 1 : 0)
              << " phi_ok=" << (phi_ok ? 1 : 0)
              << " phi_solver_rel_res=" << pipeline.em.phi_solver_rel_residual
              << " thermal_dt=" << local_cli.thermal_dt << " thermal_steps=" << local_cli.thermal_steps
              << " thermal_state_recording=" << (local_cli.thermal_state_recording ? 1 : 0)
              << " mean_delta_T=" << pipeline.thermal.mean_delta_t << " max_delta_T=" << pipeline.thermal.max_delta_t
              << " max_T=" << max_temperature << " boundary_compliance=" << boundary_compliance
              << " thermal_max_rel_err=" << pipeline.thermal.max_per_particle_rel_err
              << " vol_weighted_delta_T=" << pipeline.thermal.vol_weighted_delta_t
              << " vol_weighted_expected_delta_T=" << pipeline.thermal.vol_weighted_expected_delta_t
              << " vol_weighted_rel_err=" << vol_weighted_rel_err
              << " energy_balance_rel_err=" << pipeline.thermal.energy_balance_rel_err
              << " closure_mismatch_vol_frac=" << pipeline.thermal.closure_mismatch_vol_fraction
              << " energy_vs_power_rel_err=" << energy_vs_power_rel_err << " energy_cap_rel_err=" << energy_cap_rel_err
              << " E_joule_J=" << pipeline.thermal.total_joule_energy_j
              << " E_thermal_J=" << pipeline.thermal.total_thermal_energy_j
              << " E_power_expected_J=" << expected_energy_j << " em_ok=" << (em_ok ? 1 : 0)
              << " thermal_ok=" << (thermal_ok ? 1 : 0) << " sigma_spatial_ok=" << (sigma_spatial_ok ? 1 : 0)
              << " sigma_coupling_ok=" << (sigma_coupling_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
