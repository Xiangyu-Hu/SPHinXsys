/**
 * @file test_3d_ophelie_french_natural_em.cpp
 * @brief French natural-convection glass: complex edge-flux EM on relaxed cylinder particles.
 *
 * Geometry defaults (Jacoutot natural / EREBUS fill):
 *   R = 0.250 m, H = 0.185 m, f = 282 kHz, 7 filament loops, sigma = 16 S/m.
 *
 * Excitation:
 *   --excitation-mode=current       fixed I_per_loop (peak), report P_glass
 *   --excitation-mode=target-power  calibrate I so P_glass ≈ target (default 50 kW)
 *
 * Self-induction A_glass is intentionally off in this first EM stage.
 * Coil turn count / radius marked as thesis-supplement placeholders when not from the journal text.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_glass_mesh_relax.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_progress.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <limits>
#include <cstring>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

enum class ExcitationMode
{
    Current,
    TargetPower
};

struct LocalCli
{
    ExcitationMode excitation_mode = ExcitationMode::Current;
    bool use_reload = true;
    bool lattice_only = false;
    std::string reload_dir = "./reload";
};

inline void applyFrenchNaturalEmDefaults(OphelieParameters &params, OphelieFrenchReducedCaseParams &french)
{
    applyFrenchReducedDefaults(params, french);

    french.glass_radius = 0.250;
    french.glass_half_height = 0.5 * 0.185;
    french.glass_center = Vecd(0.0, 0.0, french.glass_half_height);
    // Stage1.5: reload/dp baseline must match relaxation spacing.
    // If Reload.xml was generated with dp=0.015, production must use dp=0.015 as well.
    french.dp = 0.015;
    french.glass_mesh_resolution = 20;

    // Thesis / EREBUS-system supplements for journal-missing coil layout.
    french.frequency_hz = 282.0e3;
    french.sigma_glass = 16.0;
    french.target_joule_power = 50.0e3;
    french.coil.num_loops = 7;
    // Stage1.5 coil geometry: EREBUS inductor diameter=570mm => coil radius ~285mm
    french.coil.loop_radius = 0.285;
    french.coil.segments_per_loop = 64;
    french.coil.use_cell_centered_loops = true;
    // Stage1.5 coil stack height total = 230mm, center-aligned with glass center.
    // glass_center[2] == glass_half_height, so half stack-height is 0.115m.
    french.coil_z_inset = 0.0;
    french.coil_z_user_set = true;
    french.coil.z_min = french.glass_center[2] - 0.115;
    french.coil.z_max = french.glass_center[2] + 0.115;
    french.auto_coil_z = true;
    french.ampere_turns = static_cast<Real>(french.coil.num_loops); // 1 A peak per loop by default
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    params.enable_phi_correction_ = true;
    params.enable_power_scaling_ = false;
    params.enable_self_induction_ = false;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::SolverLocal;
    params.phi_solver_kind_ = OpheliePhiSolverKind::GMRES;
}

inline LevelSetShape &defineGlassLevelSet(SolidBody &body)
{
#if SPHINXSYS_USE_SYCL
    return body.defineBodyLevelSetShape(par_ck, 2.0).correctLevelSetSign().cleanLevelSet();
#else
    return body.defineBodyLevelSetShape(2.0).correctLevelSetSign().cleanLevelSet();
#endif
}

inline bool parseLocalCli(int ac, char *av[], LocalCli &local_cli, StdVec<std::string> &remaining)
{
    remaining.emplace_back(av[0]);
    for (int i = 1; i < ac; ++i)
    {
        const char *arg = av[i];
        if (std::strncmp(arg, "--excitation-mode=", 18) == 0)
        {
            const std::string mode(arg + 18);
            if (mode == "current" || mode == "fixed-current" || mode == "I")
            {
                local_cli.excitation_mode = ExcitationMode::Current;
            }
            else if (mode == "target-power" || mode == "power" || mode == "P")
            {
                local_cli.excitation_mode = ExcitationMode::TargetPower;
            }
            else
            {
                std::cerr << "Unknown --excitation-mode=" << mode
                          << " (expected current|target-power)\n";
                return false;
            }
            continue;
        }
        if (std::strcmp(arg, "--lattice") == 0 || std::strcmp(arg, "--skip-reload") == 0)
        {
            local_cli.use_reload = false;
            local_cli.lattice_only = true;
            continue;
        }
        if (std::strncmp(arg, "--reload-dir=", 13) == 0)
        {
            local_cli.reload_dir = std::string(arg + 13);
            local_cli.use_reload = true;
            continue;
        }
        if (std::strcmp(arg, "--reload=1") == 0 || std::strcmp(arg, "--reload") == 0)
        {
            local_cli.use_reload = true;
            continue;
        }
        if (std::strcmp(arg, "--reload=0") == 0)
        {
            local_cli.use_reload = false;
            local_cli.lattice_only = true;
            continue;
        }
        remaining.emplace_back(arg);
    }
    return true;
}

inline void printNaturalEmSummary(const OphelieFrenchReducedCaseParams &french, ExcitationMode mode)
{
    const Real i_peak = french.coil.current_per_loop;
    const Real i_rms = i_peak / std::sqrt(Real(2));
    std::cout << "[ophelie] french_natural_em"
              << " R=" << french.glass_radius << " H=" << frenchReducedGlassHeight(french)
              << " dp=" << french.dp << " f=" << french.frequency_hz << " Hz"
              << " sigma=" << french.sigma_glass << " S/m"
              << " N_turn=" << french.coil.num_loops << " coil_R=" << french.coil.loop_radius
              << " segments=" << french.coil.segments_per_loop
              << " I_peak_per_loop=" << i_peak << " A I_rms_per_loop=" << i_rms << " A"
              << " ampere_turns=" << french.ampere_turns
              << " excitation=" << (mode == ExcitationMode::Current ? "current" : "target-power")
              << " target_P=" << french.target_joule_power << " W"
              << " edge_flux_complex=1 self_induction=0"
              << "\n[ophelie] note: N_turn/coil_R are EREBUS-system supplements when not stated in the journal text"
              << std::endl;
    logFrenchReducedCoilGeometry(french);
}

} // namespace

int main(int ac, char *av[])
{
    logOphelieRunContext();

    LocalCli local_cli;
    StdVec<std::string> after_local;
    if (!parseLocalCli(ac, av, local_cli, after_local))
    {
        return 1;
    }

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchNaturalEmDefaults(params, french);

    StdVec<char *> after_local_av;
    after_local_av.reserve(after_local.size());
    for (auto &argument : after_local)
    {
        after_local_av.push_back(const_cast<char *>(argument.c_str()));
    }

    const StdVec<std::string> french_filtered =
        filterFrenchReducedCommandLine(static_cast<int>(after_local_av.size()), after_local_av.data(), french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);

    OphelieTestCliOptions cli_options;
    StdVec<char *> french_av;
    french_av.reserve(french_filtered.size());
    for (auto &argument : french_filtered)
    {
        french_av.push_back(const_cast<char *>(argument.c_str()));
    }
    const StdVec<std::string> filtered_arguments = filterOphelieTestCommandLine(
        static_cast<int>(french_av.size()), french_av.data(), params, cli_options);

    // Force production EM path for this case.
    params.enable_phi_correction_ = true;
    params.enable_power_scaling_ = false;
    params.enable_self_induction_ = false;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    cli_options.ophelie_current_form_user_set = true;
    finalizeOphelieCurrentFormConfiguration(params, cli_options);

    french.sigma_glass = params.sigma_glass_;
    french.frequency_hz = params.frequency_;
    french.target_joule_power = params.target_joule_power_;
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);

    if (local_cli.excitation_mode == ExcitationMode::TargetPower)
    {
        OphelieFrenchLiteratureProfile literature_profile;
        literature_profile.calibrate_coil_current = true;
        applyFrenchLiteratureMode(params, cli_options, literature_profile);
        params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
        params.edge_flux_complex_ = true;
        finalizeOphelieCurrentFormConfiguration(params, cli_options);
    }

    logOphelieFinalParams(params, cli_options);
    printNaturalEmSummary(french, local_cli.excitation_mode);

    if (!cli_options.reload_dir.empty())
    {
        local_cli.reload_dir = cli_options.reload_dir;
        local_cli.use_reload = true;
    }

    const BoundingBoxd system_bounds = frenchReducedDomainBounds(french, 3.0 * french.dp);
    SPHSystem sph_system(system_bounds, french.dp);

    StdVec<char *> filtered_argv;
    filtered_argv.reserve(filtered_arguments.size());
    for (auto &argument : filtered_arguments)
    {
        filtered_argv.push_back(const_cast<char *>(argument.c_str()));
    }
    sph_system.handleCommandlineOptions(static_cast<int>(filtered_argv.size()), filtered_argv.data());

    if (local_cli.use_reload && !local_cli.lattice_only)
    {
        sph_system.setReloadParticles(true);
        IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);
        std::cout << "[ophelie] reload folder: " << IO::getEnvironment().ReloadFolder() << std::endl;
    }
    else
    {
        sph_system.setReloadParticles(false);
        std::cout << "[ophelie] lattice particles (no reload); for production use relaxed Reload.xml\n";
    }

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>(
                             "GlassBody", french.glass_center, french.glass_radius, french.glass_half_height,
                             french.glass_mesh_resolution));
    glass_body.defineAdaptation<SPHAdaptation>(1.0, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineGlassLevelSet(glass_body);

    const char *particle_source = "lattice";
    if (sph_system.ReloadParticles())
    {
        glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
        particle_source = "reload";
    }
    else
    {
        glass_body.generateParticles<BaseParticles, Lattice>();
    }

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass_fields(glass_body, glass_names);
    (void)register_glass_fields;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_glass_sigma(glass_body, glass_names,
                                                                                     params.sigma_glass_);
    syncGlassElectromagneticFieldsToDevice(glass_body.getBaseParticles(), glass_names);
    assign_glass_sigma.exec();

    BaseParticles &glass_particles = glass_body.getBaseParticles();
    const size_t n_glass = glass_particles.TotalRealParticles();

    // Stage1.5 closeout: reload/dp consistency audit (must gate).
    if (sph_system.ReloadParticles() && !local_cli.lattice_only)
    {
        syncVariableToHost<Real>(glass_particles, "VolumetricMeasure");
        const Real *vol = glass_particles.getVariableDataByName<Real>("VolumetricMeasure");

        double sum_vol = 0.0;
        double min_vol = std::numeric_limits<double>::infinity();
        double max_vol = 0.0;
        for (size_t i = 0; i < n_glass; ++i)
        {
            sum_vol += static_cast<double>(vol[i]);
            min_vol = std::min(min_vol, static_cast<double>(vol[i]));
            max_vol = std::max(max_vol, static_cast<double>(vol[i]));
        }

        const double mean_vol = sum_vol / static_cast<double>(n_glass);
        const double dp_eff = std::cbrt(mean_vol);

        const double h_exact = static_cast<double>(frenchReducedGlassHeight(french)); // = glass_half_height*2
        const double v_exact = static_cast<double>(Pi) * french.glass_radius * french.glass_radius * h_exact;
        const double vol_rel_error = std::abs(sum_vol - v_exact) / (v_exact + TinyReal);
        const double dp_rel_error = std::abs(dp_eff - french.dp) / (french.dp + TinyReal);

        std::cout << "[ophelie][stage1.5] reload audit: declared_dp=" << french.dp
                  << " dp_eff_from_meanVol=" << dp_eff << " meanVol=" << mean_vol << " minVol=" << min_vol
                  << " maxVol=" << max_vol << " sumVol=" << sum_vol << " V_exact=" << v_exact
                  << " volume_rel_error=" << vol_rel_error << " dp_rel_error=" << dp_rel_error << std::endl;

        // Neighbor-count stats for audit documentation.
        // We set phi/A fields to zeros to avoid nonfinite from uninitialized coils.
        syncVariableToHost<Real>(glass_particles, glass_names.phi_real);
        Real *phi_real = glass_particles.getVariableDataByName<Real>(glass_names.phi_real);
        syncVariableToHost<Real>(glass_particles, glass_names.phi_imag);
        Real *phi_imag = glass_particles.getVariableDataByName<Real>(glass_names.phi_imag);
        syncVariableToHost<Vecd>(glass_particles, glass_names.a_coil_real);
        Vecd *a_coil_real = glass_particles.getVariableDataByName<Vecd>(glass_names.a_coil_real);
        syncVariableToHost<Vecd>(glass_particles, glass_names.a_coil_imag);
        Vecd *a_coil_imag = glass_particles.getVariableDataByName<Vecd>(glass_names.a_coil_imag);
        for (size_t i = 0; i < n_glass; ++i)
        {
            phi_real[i] = Real(0);
            phi_imag[i] = Real(0);
            a_coil_real[i] = Vecd::Zero();
            a_coil_imag[i] = Vecd::Zero();
        }
        syncVariableToDevice<Real>(glass_particles, glass_names.phi_real);
        syncVariableToDevice<Real>(glass_particles, glass_names.phi_imag);
        syncVariableToDevice<Vecd>(glass_particles, glass_names.a_coil_real);
        syncVariableToDevice<Vecd>(glass_particles, glass_names.a_coil_imag);

        const OphelieEdgeFluxComponent imag_component = makeOphelieEdgeFluxImagComponent(glass_names, params);
        (void)evaluateOphelieEdgeFluxQAntisymmetryForComponent<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, imag_component, params);

        syncVariableToHost<Real>(glass_particles, glass_names.edge_q_neighbor_count);
        const Real *nc = glass_particles.getVariableDataByName<Real>(glass_names.edge_q_neighbor_count);
        std::vector<Real> nc_vec(n_glass);
        double nc_sum = 0.0;
        Real nc_min = nc[0];
        Real nc_max = nc[0];
        for (size_t i = 0; i < n_glass; ++i)
        {
            nc_vec[i] = nc[i];
            nc_sum += static_cast<double>(nc[i]);
            nc_min = std::min(nc_min, nc_vec[i]);
            nc_max = std::max(nc_max, nc_vec[i]);
        }
        std::sort(nc_vec.begin(), nc_vec.end());
        const Real nc_p50 = nc_vec[n_glass / 2];
        const double nc_mean = nc_sum / static_cast<double>(n_glass);

        std::cout << "[ophelie][stage1.5] neighbor-count: p50=" << nc_p50 << " mean=" << nc_mean
                  << " min=" << nc_min << " max=" << nc_max << std::endl;

        // Gate thresholds from Stage1.5 plan.
        const bool reload_ok = vol_rel_error < 0.03 && dp_rel_error < 0.05;
        if (!reload_ok)
        {
            std::cerr << "[ophelie][stage1.5] reload dp consistency FAIL: stop EM until relax/reload matches."
                      << std::endl;
            return 1;
        }
    }

    OphelieFrenchEmSolveResult em_result =
        runFrenchReducedEmPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params, french);

    if (local_cli.excitation_mode == ExcitationMode::TargetPower)
    {
        calibrateFrenchCoilCurrentToTargetPower(french, params, em_result.joule_power_raw);
        em_result =
            runFrenchReducedEmPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params, french);
    }

    const Real p_report =
        em_result.joule_power_recon_edge > TinyReal ? em_result.joule_power_recon_edge : em_result.joule_power_raw;
    const Real i_peak = french.coil.current_per_loop;
    const Real i_rms = i_peak / std::sqrt(Real(2));

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(glass_body, glass_names.a_src_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.a_src_imag);
    write_states.addToWrite<Vecd>(glass_body, glass_names.e_edge_recon_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.e_edge_recon_imag);
    write_states.addToWrite<Vecd>(glass_body, glass_names.j_edge_recon_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.j_edge_recon_imag);
    write_states.addToWrite<Real>(glass_body, glass_names.joule_heat);
    write_states.addToWrite<Real>(glass_body, glass_names.sigma);
    write_states.addToWrite<Real>(glass_body, glass_names.phi_imag);
    write_states.addToWrite<Real>(glass_body, glass_names.phi_real);
    write_states.addToWrite<Real>(glass_body, glass_names.edge_flux_residual_imag);
    write_states.addToWrite<Real>(glass_body, glass_names.edge_flux_residual_real);
    glass_body.setNewlyUpdated();
    write_states.writeToFile(0);

    const bool fields_ok = n_glass > 0 && p_report > TinyReal && std::isfinite(p_report) &&
                           std::isfinite(em_result.phi_solver_rel_residual);
    const bool power_ok = local_cli.excitation_mode == ExcitationMode::Current ||
                          std::abs(p_report - params.target_joule_power_) <=
                              0.05 * params.target_joule_power_;
    const bool passed = fields_ok && power_ok;

    std::cout << "test_3d_ophelie_french_natural_em source=" << particle_source << " n_glass=" << n_glass
              << " phi_solver_rel_res=" << em_result.phi_solver_rel_residual
              << " phi_eq_res_vol=" << em_result.phi_eq_res_vol << " P_raw=" << em_result.joule_power_raw
              << " P_recon_edge=" << em_result.joule_power_recon_edge << " P_graph_edge="
              << em_result.joule_power_graph_edge << " I_peak=" << i_peak << " I_rms=" << i_rms;
    if (em_result.edge_flux_report_valid)
    {
        std::cout << " edge_res_red=" << em_result.edge_flux_report.edge_res_red_l2;
    }
    std::cout << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
