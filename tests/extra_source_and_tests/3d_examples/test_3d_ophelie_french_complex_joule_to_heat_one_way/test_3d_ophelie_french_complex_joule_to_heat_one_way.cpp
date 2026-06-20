/**
 * @file test_3d_ophelie_french_complex_joule_to_heat_one_way.cpp
 * @brief French reload: complex edge-flux EM → JouleHeat → one-way thermal (no feedback).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_thermal_material.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"
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

    OphelieThermalDiffusionOneWayOptions thermal_options;
    thermal_options.enable_diffusion = local_cli.thermal_diffusion;
    thermal_options.enable_cold_wall_dirichlet = local_cli.thermal_diffusion;
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
    Real boundary_compliance = 1.0;
    Real max_temperature = material.t_initial;
    if (local_cli.thermal_diffusion)
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
    const bool em_ok = n > 0 && std::isfinite(pipeline.joule_power_w) && pipeline.joule_power_w > TinyReal &&
                       pipeline.phi_eq_res_vol < Real(0.01);
    const Real vol_weighted_rel_err =
        std::abs(pipeline.thermal.vol_weighted_delta_t - pipeline.thermal.vol_weighted_expected_delta_t) /
        (pipeline.thermal.vol_weighted_expected_delta_t + TinyReal);
    const Real expected_energy_j = pipeline.joule_power_w * local_cli.thermal_dt * static_cast<Real>(local_cli.thermal_steps);
    const Real energy_vs_power_rel_err =
        std::abs(pipeline.thermal.total_thermal_energy_j - expected_energy_j) / (expected_energy_j + TinyReal);
    const Real energy_cap_rel_err =
        (pipeline.thermal.total_thermal_energy_j - pipeline.thermal.total_joule_energy_j) /
        (pipeline.thermal.total_joule_energy_j + TinyReal);

    bool thermal_ok = false;
    if (local_cli.thermal_diffusion)
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
    const bool passed = em_ok && thermal_ok;

    std::cout << "test_3d_ophelie_french_complex_joule_to_heat_one_way particles=reload n=" << n
              << " dp=" << french.dp << " edge_flux_complex=1 thermal_diffusion=" << (local_cli.thermal_diffusion ? 1 : 0)
              << " thermal_material=" << ophelieThermalMaterialPresetName(local_cli.material_preset)
              << " T0=" << material.t_initial << " rho=" << material.rho << " cp=" << material.cp << " k=" << material.k
              << " P_joule_W=" << pipeline.joule_power_w << " phi_eq_res_vol=" << pipeline.phi_eq_res_vol
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
              << " thermal_ok=" << (thermal_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
