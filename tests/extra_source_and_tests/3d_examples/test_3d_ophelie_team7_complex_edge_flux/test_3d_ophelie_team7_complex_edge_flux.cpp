/**
 * @file test_3d_ophelie_team7_complex_edge_flux.cpp
 * @brief TEAM7 validation with volume-racetrack coil source + complex edge-flux (L1/L2/L3).
 *
 * Requires native STL reload (CoilSourceBody + PlateBody). No air body in EM.
 *
 * Levels (--team7-level=):
 *   coil-only  L1: coil Biot Bz vs reference (geometry/units/NI audit)
 *   one-way    L2: complex edge-flux + A_ind one-way (coil Bz gate; B_total logged for calibration)
 *   picard     L3: self-induction Picard + B_total probe
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_team7_coil_path_source.h"
#include "electromagnetic_ophelie_team7_native_geometry.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"
#include "electromagnetic_ophelie_team7_edge_recon_boundary.h"
#include "electromagnetic_ophelie_team7_boundary_consistency.h"
#include "electromagnetic_ophelie_team7_l1_source_audit.h"
#include "electromagnetic_ophelie_team7_probe.h"
#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_edge_flux_operator_audit.h"
#include "electromagnetic_ophelie_aind_lenz_audit.h"
#include "electromagnetic_ophelie_team7_validation.h"
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

constexpr const char *kTeam7DepthBinZMm = "Team7DepthBinZ_mm";

struct Team7ComplexLocalCli
{
    Team7ValidationLevel level = Team7ValidationLevel::OneWay;
    bool require_reload = true;
    std::string reload_dir;
    std::string reference_dir;
    /** Volume-racetrack vs TEAM7 filament reference (~peak Bz ratio); override with --team7-coil-source-scale=. */
    Real coil_source_scale = 0.754;
    /** 0=off; positive=factor; auto via --team7-ind-j-post-scale=auto scales J so Bind/Bcoil~1 (diagnostic only). */
    Real ind_j_post_scale = 0.0;
    bool ind_j_post_scale_auto = false;
    /** Probe-layer TEAM7 wt=90° mapping (default minus-imag: B_phase90 = −B_imag). */
    Team7Phase90Convention phase90_convention = Team7Phase90Convention::MinusImag;
    Team7ValidationMode validation_mode = Team7ValidationMode::Smoke;
    std::string output_tag;
    Real phase90_rms_abs_threshold_mT = 0.5;
    Real bz_rms_threshold_coil_only = 0.70;
    Real bz_rms_threshold_one_way = 0.40;
    Real bz_rms_threshold_picard = 0.35;
    Real team7_frequency_hz = 50.0;
    Team7BzProbeLine probe_line = Team7BzProbeLine::A1B1;
    bool team7_omega_sweep = false;
    bool team7_normalization_sweep = false;
    bool team7_operator_audit = false;
    bool team7_aind_lenz_audit = false;
    bool team7_boundary_normal_audit = false;
    bool team7_boundary_consistency_audit = false;
    bool team7_l1_source_audit = false;
};

inline bool reloadXmlExists(const std::string &folder)
{
    return fs::exists(fs::path(folder) / "Reload.xml");
}

inline std::string resolveDefaultTeam7ReloadFolder()
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
    return "./reload";
}

inline std::string resolveDefaultTeam7ReferenceDir()
{
    const StdVec<std::string> candidates = {
        "./reference_data/team7",
        "../reference_data/team7",
        "../../reference_data/team7",
        "../../../reference_data/team7",
        "../../../../../tests/extra_source_and_tests/3d_examples/reference_data/team7",
        "tests/extra_source_and_tests/3d_examples/reference_data/team7",
    };
    for (const std::string &candidate : candidates)
    {
        if (fs::exists(candidate + "/TEAM7_Bz_A1_B1_reference_mT.csv"))
        {
            return candidate;
        }
    }
    return candidates.front();
}

inline void applyTeam7ComplexLocalCli(int ac, char *av[], Team7ComplexLocalCli &local_cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            local_cli.require_reload = true;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            local_cli.reload_dir = std::string(av[i] + 13);
            local_cli.require_reload = true;
        }
        else if (std::strncmp(av[i], "--team7-reference-dir=", 22) == 0)
        {
            local_cli.reference_dir = std::string(av[i] + 22);
        }
        else if (std::strncmp(av[i], "--team7-bz-rms-threshold=", 25) == 0)
        {
            const Real threshold = static_cast<Real>(std::atof(av[i] + 25));
            local_cli.bz_rms_threshold_one_way = threshold;
            local_cli.bz_rms_threshold_picard = threshold;
            local_cli.bz_rms_threshold_coil_only = threshold;
        }
        else if (std::strncmp(av[i], "--team7-coil-source-scale=", 26) == 0)
        {
            local_cli.coil_source_scale = static_cast<Real>(std::atof(av[i] + 26));
        }
        else if (std::strncmp(av[i], "--team7-ind-j-post-scale=", 25) == 0)
        {
            const char *value = av[i] + 25;
            if (std::strcmp(value, "auto") == 0)
            {
                local_cli.ind_j_post_scale_auto = true;
            }
            else
            {
                local_cli.ind_j_post_scale = static_cast<Real>(std::atof(value));
            }
        }
        else if (std::strcmp(av[i], "--team7-probe-phase90-neg-ind=1") == 0 ||
                 std::strcmp(av[i], "--team7-probe-phase90-neg-ind") == 0)
        {
            local_cli.phase90_convention = Team7Phase90Convention::MinusImag;
        }
        else if (std::strncmp(av[i], "--team7-phase90-convention=", 27) == 0)
        {
            local_cli.phase90_convention = parseTeam7Phase90Convention(std::string(av[i] + 27));
        }
        else if (std::strncmp(av[i], "--team7-validation-mode=", 24) == 0)
        {
            local_cli.validation_mode = parseTeam7ValidationMode(std::string(av[i] + 24));
        }
        else if (std::strncmp(av[i], "--team7-output-tag=", 19) == 0)
        {
            local_cli.output_tag = std::string(av[i] + 19);
        }
        else if (std::strncmp(av[i], "--team7-phase90-rms-abs-threshold-mT=", 37) == 0)
        {
            local_cli.phase90_rms_abs_threshold_mT = static_cast<Real>(std::atof(av[i] + 37));
        }
        else if (std::strncmp(av[i], "--team7-frequency=", 18) == 0)
        {
            local_cli.team7_frequency_hz = static_cast<Real>(std::atof(av[i] + 18));
        }
        else if (std::strncmp(av[i], "--team7-probe-line=", 20) == 0)
        {
            local_cli.probe_line = parseTeam7BzProbeLine(std::string(av[i] + 20));
        }
        else if (std::strcmp(av[i], "--team7-omega-sweep=1") == 0 || std::strcmp(av[i], "--team7-omega-sweep") == 0)
        {
            local_cli.team7_omega_sweep = true;
        }
        else if (std::strcmp(av[i], "--team7-normalization-sweep=1") == 0 ||
                 std::strcmp(av[i], "--team7-normalization-sweep") == 0)
        {
            local_cli.team7_normalization_sweep = true;
        }
        else if (std::strcmp(av[i], "--team7-operator-audit=1") == 0 ||
                 std::strcmp(av[i], "--team7-operator-audit") == 0)
        {
            local_cli.team7_operator_audit = true;
        }
        else if (std::strcmp(av[i], "--team7-aind-lenz-audit=1") == 0 ||
                 std::strcmp(av[i], "--team7-aind-lenz-audit") == 0)
        {
            local_cli.team7_aind_lenz_audit = true;
        }
        else if (std::strcmp(av[i], "--team7-boundary-normal-audit=1") == 0 ||
                 std::strcmp(av[i], "--team7-boundary-normal-audit") == 0)
        {
            local_cli.team7_boundary_normal_audit = true;
        }
        else if (std::strcmp(av[i], "--team7-boundary-consistency-audit=1") == 0 ||
                 std::strcmp(av[i], "--team7-boundary-consistency-audit") == 0)
        {
            local_cli.team7_boundary_consistency_audit = true;
        }
        else if (std::strcmp(av[i], "--team7-l1-source-audit=1") == 0 ||
                 std::strcmp(av[i], "--team7-l1-source-audit") == 0)
        {
            local_cli.team7_l1_source_audit = true;
        }
    }
    (void)parseTeam7ValidationLevelCli(ac, av, local_cli.level);
    if (local_cli.level == Team7ValidationLevel::CoilOnly ||
        local_cli.level == Team7ValidationLevel::OneWay)
    {
        local_cli.team7_l1_source_audit = true;
    }
    if (local_cli.reload_dir.empty())
    {
        local_cli.reload_dir = resolveDefaultTeam7ReloadFolder();
    }
    if (local_cli.reference_dir.empty())
    {
        local_cli.reference_dir = resolveDefaultTeam7ReferenceDir();
    }
}

#if SPHINXSYS_USE_SYCL
inline LevelSetShape &defineTeam7SolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
}
#else
inline LevelSetShape &defineTeam7SolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
}
#endif

} // namespace

int main(int ac, char *av[])
{
    logOphelieRunContext();

    Team7ComplexLocalCli local_cli;
    applyTeam7ComplexLocalCli(ac, av, local_cli);

    OphelieParameters params;
    params.sigma_glass_ = 3.54e7;
    params.frequency_ = local_cli.team7_frequency_hz;
    params.target_joule_power_ = 50.0e3;
    params.enable_power_scaling_ = false;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_phi_correction_ = local_cli.level != Team7ValidationLevel::CoilOnly;
    params.enable_self_induction_ = local_cli.level == Team7ValidationLevel::Picard;
    // L2 one-way: phi from A_coil only, then A_ind; feedback re-solve with A_total is optional via CLI.
    params.aind_one_way_feedback_resolve_ = false;

    OphelieTestCliOptions cli_options;
    const StdVec<std::string> filtered_arguments = filterOphelieTestCommandLine(ac, av, params, cli_options);
    cli_options.native_stl = true;
    cli_options.no_power_scaling = true;
    if (!cli_options.coil_source_model_user_set)
    {
        cli_options.coil_source_model = OphelieCoilSourceModel::VolumeRacetrack;
    }
    finalizeOphelieCurrentFormConfiguration(params, cli_options);

    StdVec<char *> filtered_argv;
    filtered_argv.reserve(filtered_arguments.size());
    for (auto &argument : filtered_arguments)
    {
        filtered_argv.push_back(const_cast<char *>(argument.c_str()));
    }
    const int filtered_ac = static_cast<int>(filtered_argv.size());
    char **filtered_av = filtered_argv.data();

    OphelieTeam7NativeMesh native_mesh;
    native_mesh.dp_ref_mm_ = cli_options.native_dp_mm;
    native_mesh.dp_coil_mm_ = cli_options.native_dp_mm;
    native_mesh.dp_plate_mm_ = cli_options.native_dp_mm;
    native_mesh.relax_air_particles_ = false;
    if (cli_options.relaxation_steps > 0)
    {
        native_mesh.relaxation_steps_ = cli_options.relaxation_steps;
    }
    if (cli_options.relaxation_log_every > 0)
    {
        native_mesh.relaxation_log_every_ = cli_options.relaxation_log_every;
    }
    native_mesh.relaxation_vtp_every_ = cli_options.relaxation_vtp_every;
    if (cli_options.team7_coil_turns > TinyReal)
    {
        native_mesh.team7_coil_turns_ = cli_options.team7_coil_turns;
    }

    const Vec3d air_lower = native_mesh.air_lower_small_;
    const Vec3d air_upper = native_mesh.air_upper_small_;
    const Real dp = native_mesh.dp_ref_mm_ * native_mesh.stl_scale_to_meter_;
    const BoundingBoxd system_bounds(native_mesh.stl_scale_to_meter_ * air_lower,
                                     native_mesh.stl_scale_to_meter_ * air_upper);

    SPHSystem sph_system(system_bounds, dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    if (!cli_options.reload_dir.empty())
    {
        local_cli.reload_dir = cli_options.reload_dir;
    }
    IO::getEnvironment().resetReloadFolder(local_cli.reload_dir, true);
    const bool use_neg_ind_phase90 = team7Phase90ConventionUsesNegInd(local_cli.phase90_convention);
    std::cout << "[ophelie] TEAM7 complex edge-flux reload=" << local_cli.reload_dir
              << " level=" << team7ValidationLevelName(local_cli.level) << " dp_mm=" << native_mesh.dp_ref_mm_
              << " f_hz=" << local_cli.team7_frequency_hz << " probe_line=" << team7BzProbeLineName(local_cli.probe_line)
              << " imag_a_sign=" << params.edge_flux_imag_a_sign_
              << " phase90_convention=" << team7Phase90ConventionName(local_cli.phase90_convention)
              << " validation_mode=" << team7ValidationModeName(local_cli.validation_mode) << std::endl;
    std::cout << "[team7] TEAM7 validation uses fixed coil excitation (2742 AT), not target-power calibration."
              << std::endl;
    if (use_neg_ind_phase90)
    {
        std::cout << "[team7] phase90 convention: B_phase90_sim = -B_ind_imag (TEAM7 wt=90 probe mapping)"
                  << std::endl;
    }
    else
    {
        std::cout << "[team7] phase90 convention: B_phase90_sim = +B_ind_imag" << std::endl;
    }
    printTeam7DiagnosticWarnings(local_cli.coil_source_scale, params.edge_flux_imag_a_sign_,
                                 local_cli.ind_j_post_scale_auto, local_cli.ind_j_post_scale, local_cli.level,
                                 params.edge_flux_normalization_mode_, params.aind_one_way_feedback_resolve_,
                                 local_cli.team7_normalization_sweep);
    sph_system.handleCommandlineOptions(filtered_ac, filtered_av);
    const bool use_particle_reload = !sph_system.RunParticleRelaxation() && sph_system.ReloadParticles();
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles() && !reloadXmlExists(local_cli.reload_dir))
    {
        std::cerr << "error: Reload.xml not found in " << local_cli.reload_dir
                  << " (run with --relax=1 to generate, or set --reload-dir=)" << std::endl;
        return 1;
    }

    const std::string coil_stl_path = team7NativeCoilStlPath();
    const std::string plate_stl_path = team7NativePlateStlPath();
    SolidBody plate_body(sph_system,
                         makeShared<OphelieTeam7NativePlateShape>("PlateBody", plate_stl_path, native_mesh.stl_scale_to_meter_));
    SolidBody coil_body(sph_system,
                        makeShared<OphelieTeam7NativeCoilShape>("CoilSourceBody", coil_stl_path, native_mesh.stl_scale_to_meter_));
    plate_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    coil_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    plate_body.defineMatterMaterial<Solid>();
    coil_body.defineMatterMaterial<Solid>();
    LevelSetShape &plate_level_set = defineTeam7SolidLevelSet(plate_body);
    LevelSetShape &coil_level_set = defineTeam7SolidLevelSet(coil_body);

    if (sph_system.RunParticleRelaxation())
    {
#if SPHINXSYS_USE_SYCL
        plate_body.generateParticles<BaseParticles, Lattice>();
        coil_body.generateParticles<BaseParticles, Lattice>();
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        relaxTeam7NativeStlBodies(sph_system, coil_body, plate_body, coil_level_set, plate_level_set, native_mesh,
                                  nullptr, nullptr);
        computeTeam7CoilPlateNormalsFromShape(sph_system, coil_body, plate_body);
        ReloadParticleIO write_reload({&coil_body, &plate_body});
        registerTeam7CoilPlateNormalsForReload(write_reload, coil_body, plate_body);
        write_reload.writeToFile(0);
        std::cout << "[team7] particle relaxation finished with SPHinXsys NormalDirection -> "
                  << IO::getEnvironment().ReloadFolder() << "/Reload.xml" << std::endl;
        if (!sph_system.ReloadParticles())
        {
            return 0;
        }
        std::cout << "[team7] continuing to EM run after --relax=1 --reload=1" << std::endl;
        reloadTeam7CoilPlateParticlesWithNormals(coil_body, plate_body);
#else
        std::cerr << "error: TEAM7 native --relax=1 requires SPHINXSYS_USE_SYCL build" << std::endl;
        return 1;
#endif
    }
    else if (use_particle_reload)
    {
        reloadTeam7CoilPlateParticlesWithNormals(coil_body, plate_body);
        std::cout << "[ophelie] loaded relaxed particles + NormalDirection from "
                  << IO::getEnvironment().ReloadFolder() << "/Reload.xml" << std::endl;
    }
    else
    {
        std::cerr << "error: TEAM7 complex edge-flux requires --reload=1 (default) or --relax=1 [--reload=1 to continue]"
                  << std::endl;
        return 1;
    }

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    if (!team7CoilPlateHaveSphinxNormals(coil_body, plate_body))
    {
        std::cerr << "[team7] warning: NormalDirection/SignedDistance missing after reload; "
                     "regenerate reload with --relax=1"
                  << std::endl;
    }
    else if (local_cli.team7_boundary_normal_audit || local_cli.team7_operator_audit)
    {
        const std::string output_tag = local_cli.output_tag.empty() ? "default" : local_cli.output_tag;
        writeTeam7BoundaryNormalAuditCsv(
            coil_body, plate_body, dp,
            team7OutputArtifactPath("team7_boundary_normal_audit_", local_cli.level, output_tag, ".csv"));
    }

    const OphelieTeam7NativeDerivedGeometry native_derived = deriveTeam7NativeGeometry(coil_body, plate_body);
    applyTeam7NativeParameters(params, native_mesh, native_derived, coil_body.getBaseParticles());
    params.frequency_ = local_cli.team7_frequency_hz;

    OphelieTeam7CoilPathSourceSpec coil_path_spec;
    coil_path_spec.turns = native_mesh.team7_coil_turns_;
    coil_path_spec.current_per_turn = native_mesh.team7_coil_current_per_turn_;
    coil_path_spec.source_scale = local_cli.coil_source_scale;
    const bool use_filament_source = cli_options.coil_source_model == OphelieCoilSourceModel::FilamentRacetrack;
    const Real filament_ds_m = cli_options.racetrack_ds_mm * 1.0e-3;
    OphelieTeam7CoilPathPrepareSummary coil_path_summary;
    OphelieTeam7FilamentPrepareSummary filament_summary;
    StdVec<OphelieCurrentMomentSample> filament_moments;
    if (use_filament_source)
    {
        filament_summary =
            prepareOphelieTeam7FilamentRacetrackCoilSource(native_derived.coil_bbox_, coil_path_spec, filament_ds_m);
        buildOphelieTeam7FilamentMomentsFromCoilPath(native_derived.coil_bbox_, coil_path_spec, filament_ds_m,
                                                     filament_moments);
        printOphelieTeam7FilamentPrepareSummary(filament_summary, coil_path_spec);
        params.coil_j0_override_ = 0.0;
    }
    else
    {
        coil_path_summary =
            prepareOphelieTeam7VolumeRacetrackCoilSource(coil_body, native_derived.coil_bbox_, coil_path_spec);
        params.coil_j0_override_ = coil_path_summary.j0;
        printOphelieTeam7CoilPathPrepareSummary(coil_path_summary, coil_path_spec);
    }
    std::cout << "[ophelie] TEAM7 coil_source_model=" << ophelieCoilSourceModelName(cli_options.coil_source_model)
              << std::endl;

    OphelieCoilFieldNames coil_names;
    OphelieGlassFieldNames plate_names;
    RegisterOphelieCoilFields register_coil(coil_body, coil_names);
    RegisterOphelieGlassFields register_plate(plate_body, plate_names);
    (void)register_coil;
    (void)register_plate;
    plate_body.getBaseParticles().registerStateVariable<Real>(kTeam7DepthBinZMm, Real(0));
    UniquePtr<Inner<>> plate_inner = makeUnique<Inner<>>(plate_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(plate_body, plate_names, params.sigma_glass_);
    assign_sigma.exec();

    bool coil_path_audit_ok = false;
    if (use_filament_source)
    {
        applyOphelieTeam7FilamentRacetrackBiotToGlass(plate_body, plate_names, filament_moments, params.mu0_,
                                                      params.softening_length_);
        coil_path_audit_ok = ophelieTeam7FilamentAmpereTurnsAuditPassed(filament_summary, coil_path_spec.source_scale);
    }
    else
    {
        StateDynamics<MainExecutionPolicy, InitializeOphelieVolumeRacetrackCoilSourceCK> init_coil_source(
            coil_body, coil_names, coil_path_summary.j0);
        init_coil_source.exec();
        const Real integrated_current =
            hostCoilIntegratedCurrentFromJSrc(coil_body.getBaseParticles(), coil_names, coil_path_summary.path_length_m);
        coil_path_audit_ok = ophelieTeam7CoilPathAmpereTurnsAuditPassed(coil_path_summary, integrated_current,
                                                                        coil_path_spec.source_scale);
        StateDynamics<MainExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot(
            plate_body, coil_body, plate_names, coil_names, params);
        compute_biot.exec();
    }

    Team7OneWayEdgeFluxSummary one_way_summary;
    Real phi_solver_rel_residual = 0.0;
    Real self_induction_j_rel_change = 0.0;
    size_t self_induction_iterations_used = 0;
    bool picard_converged = false;

    StdVec<Team7BzProbePoint> reference_probes;
    const bool reference_ok = loadTeam7BzReference(local_cli.reference_dir, local_cli.probe_line,
                                                   local_cli.team7_frequency_hz, 1.0e-3, reference_probes);
    const Real plate_z_top_m = native_derived.plate_bbox_.upper_[2];
    const Real plate_z_mid_m =
        Real(0.5) * (native_derived.plate_bbox_.lower_[2] + native_derived.plate_bbox_.upper_[2]);
    const Real plate_z_skin_min_m = plate_z_top_m - Real(2) * dp;
    const StdVec<OphelieCurrentMomentSample> *filament_moments_ptr =
        use_filament_source && !filament_moments.empty() ? &filament_moments : nullptr;

    Team7EdgeFluxOmegaScalingRecord omega_scaling_record;
    if (local_cli.level == Team7ValidationLevel::CoilOnly)
    {
        std::cout << "[ophelie] TEAM7 L1 coil-only: skip phi / edge-flux / A_ind" << std::endl;
    }
    else
    {
        if (local_cli.team7_omega_sweep)
        {
            const StdVec<Real> sweep_frequencies_hz = {Real(25), Real(50), Real(100), Real(200)};
            std::cout << "[team7] edge-flux omega scaling sweep (diagnostic): n_freq=" << sweep_frequencies_hz.size()
                      << std::endl;
            const StdVec<Team7EdgeFluxOmegaScalingRecord> sweep_records =
                runTeam7EdgeFluxOmegaScalingSweep<MainExecutionPolicy>(
                    plate_body, *plate_inner, plate_names, params, dp, sweep_frequencies_hz, params.sigma_glass_,
                    filament_moments, use_filament_source, coil_body, coil_names, coil_path_summary);
            const std::string sweep_tag =
                local_cli.output_tag.empty()
                    ? makeTeam7AutoOutputTag(local_cli.level, params.frequency_, params.edge_flux_imag_a_sign_,
                                             local_cli.phase90_convention, local_cli.coil_source_scale, false, 0.0,
                                             cli_options.coil_source_model)
                    : local_cli.output_tag;
            writeTeam7EdgeFluxOmegaScalingCsv(
                team7OutputArtifactPath("team7_omega_scaling_", local_cli.level, sweep_tag, ".csv"), sweep_records);
        }
        params.frequency_ = local_cli.team7_frequency_hz;
        if (use_filament_source)
        {
            applyOphelieTeam7FilamentRacetrackBiotToGlass(plate_body, plate_names, filament_moments, params.mu0_,
                                                          params.softening_length_);
        }
        else
        {
            StateDynamics<MainExecutionPolicy, InitializeOphelieVolumeRacetrackCoilSourceCK> init_coil_source(
                coil_body, coil_names, coil_path_summary.j0);
            init_coil_source.exec();
            StateDynamics<MainExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot(
                plate_body, coil_body, plate_names, coil_names, params);
            compute_biot.exec();
        }
    }
    if (local_cli.team7_normalization_sweep && local_cli.level != Team7ValidationLevel::CoilOnly)
    {
        const StdVec<Real> safe_rhs_l2_values = {Real(1.0e4), Real(1.0e5), Real(1.0e6), Real(1.0e12)};
        std::cout << "[team7] edge-flux normalization sweep (field-scale-restore): n_safe_rhs="
                  << safe_rhs_l2_values.size() << std::endl;
        const StdVec<Team7NormalizationSweepRecord> norm_records =
            runTeam7EdgeFluxNormalizationSweep<MainExecutionPolicy>(
                plate_body, *plate_inner, plate_names, params, dp, safe_rhs_l2_values, reference_probes, reference_ok,
                local_cli.phase90_convention, plate_z_skin_min_m, plate_z_top_m, plate_z_mid_m, filament_moments,
                use_filament_source, coil_body, coil_names, coil_path_summary);
        const std::string norm_tag =
            local_cli.output_tag.empty()
                ? makeTeam7AutoOutputTag(local_cli.level, params.frequency_, params.edge_flux_imag_a_sign_,
                                         local_cli.phase90_convention, local_cli.coil_source_scale, false, 0.0,
                                         cli_options.coil_source_model)
                : local_cli.output_tag;
        writeTeam7NormalizationSweepCsv(
            team7OutputArtifactPath("team7_normalization_sweep_", local_cli.level, norm_tag, ".csv"), norm_records);
        if (use_filament_source)
        {
            applyOphelieTeam7FilamentRacetrackBiotToGlass(plate_body, plate_names, filament_moments, params.mu0_,
                                                          params.softening_length_);
        }
        else
        {
            StateDynamics<MainExecutionPolicy, InitializeOphelieVolumeRacetrackCoilSourceCK> init_coil_source(
                coil_body, coil_names, coil_path_summary.j0);
            init_coil_source.exec();
            StateDynamics<MainExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot(
                plate_body, coil_body, plate_names, coil_names, params);
            compute_biot.exec();
        }
    }
    if (local_cli.level == Team7ValidationLevel::OneWay)
    {
        one_way_summary =
            runTeam7ComplexEdgeFluxOneWay<MainExecutionPolicy>(plate_body, *plate_inner, plate_names, params, dp);
        omega_scaling_record =
            makeTeam7EdgeFluxOmegaScalingRecord(local_cli.team7_frequency_hz, params.sigma_glass_, one_way_summary);
        printTeam7EdgeFluxOmegaScalingReport({omega_scaling_record}, local_cli.team7_frequency_hz);
    }

    Team7BzCompareMetrics coil_metrics;
    Team7BzCompareMetrics coil_span_metrics;
    Team7BzCompareMetrics total_metrics;
    Team7BzCompareMetrics total_phase90_metrics;
    Team7BzCompareMetrics total_magnitude_metrics;
    Team7BzCompareMetrics ind_imag_particle_metrics;
    Team7BzCompareMetrics ind_imag_skin_metrics;
    StdVec<Team7ProbeXBandMetrics> x_band_metrics;
    coil_metrics.passed = true;
    coil_span_metrics.passed = true;
    total_metrics.passed = true;
    total_phase90_metrics.passed = true;
    total_magnitude_metrics.passed = true;
    ind_imag_particle_metrics.passed = true;
    ind_imag_skin_metrics.passed = true;

    const std::string &plate_j_real = getOphelieAIndJRealFieldName(plate_names, params);
    const std::string &plate_j_imag = getOphelieAIndJImagFieldName(plate_names, params);
    const size_t n_plate_pre_probe = plate_body.getBaseParticles().TotalRealParticles();
    Real applied_ind_j_post_scale = 0.0;
    if (local_cli.level == Team7ValidationLevel::OneWay)
    {
        const Team7ImagChainPowerAudit power_audit =
            auditTeam7ImagChainPowerOnPlate(plate_body.getBaseParticles(), plate_names, plate_j_imag, n_plate_pre_probe);
        printTeam7ImagChainPowerAudit(power_audit);
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list_jaudit(plate_body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation_jaudit(*plate_inner);
        update_cell_linked_list_jaudit.exec();
        update_inner_relation_jaudit.exec();
        InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi_jaudit(
            *plate_inner, plate_names);
        compute_grad_phi_jaudit.exec();
        for (size_t ip = 0; ip < reference_probes.size(); ++ip)
        {
            if (std::abs(reference_probes[ip].x_mm - Real(162)) < Real(0.5))
            {
                const Team7ImagJEmVsEdgeAudit j_audit = auditTeam7ImagJEmVsEdgeAtProbe(
                    plate_body.getBaseParticles(), plate_names, params, plate_j_imag, reference_probes[ip]);
                printTeam7ImagJEmVsEdgeAudit(j_audit);
                break;
            }
        }
        Real j_post_scale = 0.0;
        if (local_cli.ind_j_post_scale_auto)
        {
            j_post_scale = Real(1) / (one_way_summary.b_ind_over_b_coil + TinyReal);
        }
        else if (local_cli.ind_j_post_scale > TinyReal)
        {
            j_post_scale = local_cli.ind_j_post_scale;
        }
        if (j_post_scale > TinyReal)
        {
            applied_ind_j_post_scale = applyTeam7InducedJPostScaleAndRefreshAInd<MainExecutionPolicy>(
                plate_body, plate_names, params, j_post_scale, n_plate_pre_probe);
            const Team7ImagChainPowerAudit power_after =
                auditTeam7ImagChainPowerOnPlate(plate_body.getBaseParticles(), plate_names, plate_j_imag,
                                                n_plate_pre_probe);
            printTeam7ImagChainPowerAudit(power_after);
        }
    }

    const std::string level_tag = team7ValidationLevelName(local_cli.level);
    std::string output_tag = local_cli.output_tag;
    if (output_tag.empty())
    {
        output_tag = makeTeam7AutoOutputTag(local_cli.level, params.frequency_, params.edge_flux_imag_a_sign_,
                                            local_cli.phase90_convention, local_cli.coil_source_scale,
                                            local_cli.ind_j_post_scale_auto, applied_ind_j_post_scale,
                                            cli_options.coil_source_model);
    }
    std::cout << "[team7] output_tag=" << output_tag << std::endl;

    if (local_cli.level == Team7ValidationLevel::Picard)
    {
        const std::string picard_csv_path =
            team7OutputArtifactPath("team7_picard_", local_cli.level, output_tag, ".csv");
        const Team7PicardRunResult picard_result = runTeam7ComplexEdgeFluxPicardWithLog<MainExecutionPolicy>(
            plate_body, *plate_inner, plate_names, params, dp, coil_body.getBaseParticles(), coil_names,
            reference_probes, reference_ok, local_cli.phase90_convention, filament_moments_ptr, plate_z_skin_min_m,
            plate_z_top_m, plate_z_mid_m, picard_csv_path);
        one_way_summary = picard_result.final_summary;
        self_induction_j_rel_change = picard_result.final_j_rel_change;
        self_induction_iterations_used = picard_result.iterations_used;
        phi_solver_rel_residual = picard_result.phi_solver_rel_residual;
        picard_converged = picard_result.picard_converged;
    }

    if (local_cli.team7_operator_audit &&
        (local_cli.level == Team7ValidationLevel::OneWay || local_cli.level == Team7ValidationLevel::Picard))
    {
        runTeam7EdgeFluxOperatorAudit<MainExecutionPolicy>(
            plate_body, *plate_inner, plate_names, params, plate_j_imag, native_derived.plate_bbox_,
            native_derived.plate_bbox_.lower_[2], native_derived.plate_bbox_.upper_[2], dp,
            team7OutputArtifactPath("team7_edge_conductance_audit_", local_cli.level, output_tag, ".csv"),
            team7OutputArtifactPath("team7_edge_partition_audit_", local_cli.level, output_tag, ".csv"));
        local_cli.team7_aind_lenz_audit = true;
    }

    if (local_cli.team7_aind_lenz_audit &&
        (local_cli.level == Team7ValidationLevel::OneWay || local_cli.level == Team7ValidationLevel::Picard))
    {
        runTeam7AIndLenzAuditAfterOneWay(
            plate_body, plate_names, plate_j_imag,
            team7OutputArtifactPath("team7_aind_lenz_audit_", local_cli.level, output_tag, ".csv"), output_tag);
    }

    if (ophelieEdgeReconBoundaryModeIsDiagnosticOnly(params.edge_recon_boundary_mode_) &&
        (local_cli.level == Team7ValidationLevel::OneWay || local_cli.level == Team7ValidationLevel::Picard))
    {
        StdVec<Team7EdgeReconBoundaryPartitionCompare> boundary_compare_records =
            runTeam7EdgeReconBoundaryDiagnostic<MainExecutionPolicy>(
                plate_body, *plate_inner, plate_names, params, native_derived.plate_bbox_,
                native_derived.plate_bbox_.lower_[2], native_derived.plate_bbox_.upper_[2], dp,
                team7OutputArtifactPath("team7_edge_recon_boundary_", local_cli.level, output_tag, ".csv"));

        if (local_cli.team7_boundary_consistency_audit || local_cli.team7_operator_audit)
        {
            runTeam7BoundaryConsistencyAudit<MainExecutionPolicy>(
                plate_body, *plate_inner, plate_names, params, native_derived.plate_bbox_,
                native_derived.plate_bbox_.lower_[2], native_derived.plate_bbox_.upper_[2], dp,
                team7OutputArtifactPath("team7_boundary_consistency_", local_cli.level, output_tag, ".csv"),
                &boundary_compare_records);
        }
    }
    else if ((local_cli.team7_boundary_consistency_audit || local_cli.team7_operator_audit ||
              ophelieEdgeReconBoundaryModeUsesProductionClosure(params.edge_recon_boundary_mode_)) &&
             (local_cli.level == Team7ValidationLevel::OneWay || local_cli.level == Team7ValidationLevel::Picard))
    {
        runTeam7BoundaryConsistencyAudit<MainExecutionPolicy>(
            plate_body, *plate_inner, plate_names, params, native_derived.plate_bbox_,
            native_derived.plate_bbox_.lower_[2], native_derived.plate_bbox_.upper_[2], dp,
            team7OutputArtifactPath("team7_boundary_consistency_", local_cli.level, output_tag, ".csv"), nullptr);
    }

    const StdVec<Team7ProbeBzDecomposition> decomp = evaluateTeam7ProbeBzDecomposition(
        coil_body.getBaseParticles(), coil_names, plate_body.getBaseParticles(), plate_names, params, plate_j_real,
        plate_j_imag, reference_probes, plate_z_skin_min_m, plate_z_top_m, plate_z_mid_m, dp, filament_moments_ptr);

    const std::string csv_path = team7OutputArtifactPath("team7_probe_", local_cli.level, output_tag, ".csv");
    writeTeam7LevelProbeCsv(csv_path, reference_probes, decomp, reference_ok);

    Real bz_threshold_coil = local_cli.bz_rms_threshold_coil_only;
    Real bz_threshold_total = local_cli.bz_rms_threshold_coil_only;
    if (local_cli.level == Team7ValidationLevel::OneWay)
    {
        bz_threshold_total = local_cli.bz_rms_threshold_one_way;
    }
    else if (local_cli.level == Team7ValidationLevel::Picard)
    {
        bz_threshold_total = local_cli.bz_rms_threshold_picard;
    }

    if (reference_ok)
    {
        const StdVec<Team7BzProbePoint> coil_probes =
            team7ProbesFromDecompositionReal(reference_probes, decomp, false);
        const StdVec<Team7BzProbePoint> total_probes = team7ProbesFromDecompositionReal(reference_probes, decomp, true);
        coil_metrics = compareTeam7BzPhase0(coil_probes, bz_threshold_coil);
        const Real coil_x_min_mm = native_derived.coil_bbox_.lower_[0] * 1000.0;
        coil_span_metrics = compareTeam7BzPhase0OverXRange(coil_probes, coil_x_min_mm, Team7ReferenceProbeLineMm::x_end_mm,
                                                             bz_threshold_coil);
        total_metrics = compareTeam7BzPhase0(total_probes, bz_threshold_total);
        StdVec<Real> phase90_sim_mT;
        StdVec<Real> phase90_ref_mT;
        team7ProbePhase90SimAndRef(reference_probes, decomp, phase90_sim_mT, phase90_ref_mT, use_neg_ind_phase90);
        total_phase90_metrics =
            compareTeam7BzAgainstReference(reference_probes, phase90_sim_mT, phase90_ref_mT, bz_threshold_total);
        StdVec<Real> magnitude_sim_mT;
        StdVec<Real> magnitude_ref_mT;
        team7ProbeComplexMagnitudeSimAndRef(reference_probes, decomp, magnitude_sim_mT, magnitude_ref_mT);
        total_magnitude_metrics =
            compareTeam7BzAgainstReference(reference_probes, magnitude_sim_mT, magnitude_ref_mT, bz_threshold_total);
        StdVec<Real> ind_imag_particle_sim_mT;
        for (size_t i = 0; i < decomp.size(); ++i)
        {
            ind_imag_particle_sim_mT.push_back(decomp[i].bz_ind_imag_particle_mT);
        }
        ind_imag_particle_metrics = compareTeam7BzAgainstReference(reference_probes, ind_imag_particle_sim_mT,
                                                                 phase90_ref_mT, bz_threshold_total);
        StdVec<Real> ind_imag_skin_sim_mT;
        for (size_t i = 0; i < decomp.size(); ++i)
        {
            ind_imag_skin_sim_mT.push_back(decomp[i].bz_ind_imag_skin_mT);
        }
        ind_imag_skin_metrics = compareTeam7BzAgainstReference(reference_probes, ind_imag_skin_sim_mT, phase90_ref_mT,
                                                               bz_threshold_total);
        if (local_cli.level == Team7ValidationLevel::OneWay)
        {
            const Team7ProbePhase90BiotVariantMetrics biot_variants =
                evaluateTeam7ProbePhase90BiotVariants(reference_probes, decomp, bz_threshold_total);
            printTeam7ProbePhase90BiotVariantReport(biot_variants);
        }
        const Team7L1SourceValidationReport l1_report = makeTeam7L1SourceValidationReport(
            cli_options.coil_source_model, local_cli.coil_source_scale, coil_probes, coil_metrics, coil_span_metrics);
        printTeam7L1SourceValidationReport(l1_report);
        if (local_cli.team7_l1_source_audit)
        {
            Team7L1SourceAuditSummary l1_audit_summary;
            const StdVec<Team7ProbeBzDecomposition> *decomp_ptr =
                local_cli.level == Team7ValidationLevel::CoilOnly ? nullptr : &decomp;
            const StdVec<Team7L1SourceAuditRecord> l1_audit_records = computeTeam7L1SourceAudit(
                reference_probes, coil_probes, decomp_ptr, reference_ok, cli_options.coil_source_model,
                local_cli.coil_source_scale, local_cli.probe_line, native_derived, native_mesh, coil_path_spec, params,
                coil_path_audit_ok, coil_metrics, filament_ds_m, l1_audit_summary);
            printTeam7L1SourceAuditReport(l1_audit_summary, l1_audit_records);
            writeTeam7L1SourceAuditCsv(
                team7OutputArtifactPath("team7_l1_source_audit_", local_cli.level, output_tag, ".csv"),
                team7OutputArtifactPath("team7_l1_source_audit_summary_", local_cli.level, output_tag, ".csv"),
                l1_audit_summary, l1_audit_records);
        }
        if (local_cli.level != Team7ValidationLevel::CoilOnly)
        {
            const Real edge_fallback_frac =
                hostOphelieEdgeFluxFallbackFraction(plate_body.getBaseParticles(), plate_names);
            std::cout << "[ophelie] TEAM7 Bz total phase0 RMS=" << total_metrics.rms_rel_error
                      << " peak_sim/ref=" << total_metrics.peak_sim_mT << "/" << total_metrics.peak_ref_mT
                      << " phase90_RMS=" << total_phase90_metrics.rms_rel_error
                      << " peak_sim/ref=" << total_phase90_metrics.peak_sim_mT << "/"
                      << total_phase90_metrics.peak_ref_mT << " |B|_RMS=" << total_magnitude_metrics.rms_rel_error
                      << " ind_imag_skin_RMS=" << ind_imag_skin_metrics.rms_rel_error
                      << " ind_imag_particle_RMS=" << ind_imag_particle_metrics.rms_rel_error
                      << " edge_fallback_frac=" << edge_fallback_frac
                      << " edge_flux_scale=" << one_way_summary.edge_flux_input_scale
                      << " J_imag_L2_vol=" << one_way_summary.j_imag_vol_norm
                      << " P_recon_W=" << one_way_summary.joule_power_recon_w
                      << " Bind/Bcoil=" << one_way_summary.b_ind_over_b_coil;
            if (applied_ind_j_post_scale > TinyReal)
            {
                std::cout << " j_post_scale=" << applied_ind_j_post_scale
                          << " phase90_RMS_after_jscale=" << total_phase90_metrics.rms_rel_error
                          << " ind_imag_skin_RMS_after_jscale=" << ind_imag_skin_metrics.rms_rel_error;
            }
            std::cout << std::endl;
        }
        writeTeam7BzCompareCsv(team7OutputArtifactPath("team7_bz_", local_cli.level, output_tag, "_coil_only.csv"),
                               coil_probes);
        if (local_cli.level != Team7ValidationLevel::CoilOnly)
        {
            writeTeam7BzCompareCsv(team7OutputArtifactPath("team7_bz_", local_cli.level, output_tag, "_total.csv"),
                                   total_probes);

            StdVec<Team7JeyProbePoint> jey_probes;
            const bool jey_reference_ok =
                loadTeam7JeyReference(local_cli.reference_dir, local_cli.team7_frequency_hz, 1.0e-3, jey_probes);
            if (!jey_reference_ok)
            {
                buildTeam7JeyProbeLineFromBzReference(reference_probes, 1.0e-3, jey_probes);
            }
            if (!jey_probes.empty())
            {
                const Real jy_kernel_h = params.softening_length_ * Real(4);
                StdVec<Team7JeyProbePoint> jey_sim_probes;
                evaluatePlateJyAtProbesSphInterp(plate_body.getBaseParticles(), plate_j_real, plate_j_imag, jy_kernel_h,
                                                 jey_probes, jey_sim_probes);
                for (size_t i = 0; i < jey_probes.size(); ++i)
                {
                    jey_probes[i].jy_sim_phase0_Am2 = jey_sim_probes[i].jy_sim_phase0_Am2;
                    jey_probes[i].jy_sim_phase90_Am2 = jey_sim_probes[i].jy_sim_phase90_Am2;
                }
                StdVec<Real> jy_phase0_sim;
                StdVec<Real> jy_phase0_ref;
                StdVec<Real> jy_phase90_sim;
                StdVec<Real> jy_phase90_ref;
                for (const Team7JeyProbePoint &probe : jey_probes)
                {
                    jy_phase0_sim.push_back(probe.jy_sim_phase0_Am2);
                    jy_phase0_ref.push_back(probe.jy_ref_phase0_Am2);
                    jy_phase90_sim.push_back(probe.jy_sim_phase90_Am2);
                    jy_phase90_ref.push_back(probe.jy_ref_phase90_Am2);
                }
                const bool jey_frequency_ref_ok =
                    jey_reference_ok && team7JeyReferenceHasFrequencyData(jey_probes, local_cli.team7_frequency_hz);
                const Team7BzCompareMetrics jy_phase0_metrics =
                    jey_frequency_ref_ok
                        ? compareTeam7ScalarAgainstReference(jy_phase0_sim, jy_phase0_ref, Real(1.0))
                        : Team7BzCompareMetrics{};
                const Team7BzCompareMetrics jy_phase90_metrics =
                    jey_frequency_ref_ok
                        ? compareTeam7ScalarAgainstReference(jy_phase90_sim, jy_phase90_ref, Real(1.0))
                        : Team7BzCompareMetrics{};
                const Team7JeySurfaceDepthDiagnostic jey_surface_depth = computeTeam7JeySurfaceDepthDiagnostic(
                    plate_body.getBaseParticles(), plate_names, plate_j_real, plate_j_imag,
                    native_derived.plate_bbox_.lower_[2], native_derived.plate_bbox_.upper_[2], dp, dp);
                printTeam7JeyValidationReport(local_cli.team7_frequency_hz, jey_reference_ok, jey_frequency_ref_ok,
                                              jy_phase0_metrics, jy_phase90_metrics, jey_surface_depth);
                writeTeam7JeyProbeCsv(team7OutputArtifactPath("team7_jey_probe_", local_cli.level, output_tag, ".csv"),
                                      jey_probes, jey_reference_ok);
                if (jey_frequency_ref_ok)
                {
                    const StdVec<Team7JvsBProbeSplitRecord> jvsb_records = computeTeam7JvsBProbeSplit(
                        jey_probes, reference_probes, decomp, use_neg_ind_phase90);
                    printTeam7JvsBProbeSplitReport(jvsb_records, Real(162));
                    writeTeam7JvsBProbeSplitCsv(
                        team7OutputArtifactPath("team7_j_vs_b_probe_split_", local_cli.level, output_tag, ".csv"),
                        jvsb_records);
                }
            }

            const Real coil_x_min_mm = native_derived.coil_bbox_.lower_[0] * 1000.0;
            const Real coil_x_max_mm = native_derived.coil_bbox_.upper_[0] * 1000.0;
            const StdVec<Team7ProbeXBandSpec> x_bands = defaultTeam7ProbePhase90XBands(coil_x_min_mm, coil_x_max_mm);
            x_band_metrics = evaluateTeam7ProbePhase90XBands(reference_probes, decomp, x_bands, bz_threshold_total);
            printTeam7ProbePhase90XBandReport(x_band_metrics);
            writeTeam7ProbePhase90XBandCsv(
                team7OutputArtifactPath("team7_probe_phase90_xbands_", local_cli.level, output_tag, ".csv"),
                x_band_metrics);

            StdVec<Real> depth_bin_center_mm;
            const StdVec<Team7PlateDepthBinStats> depth_bins = computeTeam7PlateDepthProfile(
                plate_body.getBaseParticles(), plate_names, plate_j_imag, native_derived.plate_bbox_.lower_[2],
                native_derived.plate_bbox_.upper_[2], dp, 1.0e-3, &depth_bin_center_mm);
            hostWriteScalarField(plate_body.getBaseParticles(), kTeam7DepthBinZMm, depth_bin_center_mm);
            printTeam7PlateDepthProfileReport(depth_bins);
            writeTeam7PlateDepthProfileCsv(
                team7OutputArtifactPath("team7_plate_depth_profile_", local_cli.level, output_tag, ".csv"), depth_bins);
        }
    }

    const size_t n_plate = plate_body.getBaseParticles().TotalRealParticles();
    const Real max_j_imag =
        local_cli.level == Team7ValidationLevel::CoilOnly
            ? 0.0
            : hostVecdFieldMax(plate_body.getBaseParticles(), plate_j_imag, n_plate);
    const Real joule_p_report = one_way_summary.p_complex_total > TinyReal ? one_way_summary.p_complex_total
                                                                           : one_way_summary.p_complex_coil_only;
    const bool em_ok =
        local_cli.level == Team7ValidationLevel::CoilOnly ||
        (max_j_imag > TinyReal && std::isfinite(max_j_imag) && std::isfinite(joule_p_report) &&
         joule_p_report > TinyReal &&
         (local_cli.level != Team7ValidationLevel::OneWay || one_way_summary.phi_eq_res_vol_imag < Real(0.01)));
    const bool picard_ok = local_cli.level != Team7ValidationLevel::Picard ||
                           self_induction_j_rel_change < params.self_induction_j_tolerance_;
    const Team7ValidationPassReport pass_report = evaluateTeam7ValidationPassReport(
        local_cli.level, local_cli.validation_mode, reference_ok, coil_path_audit_ok, em_ok, picard_ok, coil_metrics,
        total_metrics, total_phase90_metrics, local_cli.coil_source_scale, params.edge_flux_imag_a_sign_,
        local_cli.ind_j_post_scale_auto, applied_ind_j_post_scale, params.edge_flux_normalization_mode_,
        params.aind_one_way_feedback_resolve_, local_cli.phase90_rms_abs_threshold_mT);
    printTeam7ValidationPassReport(pass_report);

    std::cout << "test_3d_ophelie_team7_complex_edge_flux level=" << level_tag << " tag=" << output_tag
              << " n_coil=" << coil_body.getBaseParticles().TotalRealParticles()
              << " n_plate=" << plate_body.getBaseParticles().TotalRealParticles()
              << " coil_path_audit=" << (coil_path_audit_ok ? 1 : 0)
              << " bz_rms_coil=" << coil_metrics.rms_rel_error
              << " bz_rms_coil_x_span=" << coil_span_metrics.rms_rel_error
              << " bz_rms_total=" << total_metrics.rms_rel_error
              << " bz_rms_total_phase90=" << total_phase90_metrics.rms_rel_error
              << " phase90_rms_abs_mT=" << total_phase90_metrics.rms_abs_error_mT
              << " bz_rms_total_mag=" << total_magnitude_metrics.rms_rel_error
              << " bz_rms_ind_imag_skin=" << ind_imag_skin_metrics.rms_rel_error
              << " smoke_passed=" << (pass_report.smoke_passed ? 1 : 0)
              << " team7_validation_passed=" << (pass_report.team7_validation_passed ? 1 : 0)
              << " diagnostic_only=" << (pass_report.diagnostic_only ? 1 : 0) << std::endl;

    if (reference_ok)
    {
        const std::string history_path =
            team7OutputArtifactPath("team7_validation_history_", local_cli.level, output_tag, ".txt");
        appendTeam7ValidationRunRecord(history_path, level_tag, output_tag, pass_report, coil_metrics, coil_span_metrics,
                                       total_metrics, total_phase90_metrics, total_magnitude_metrics,
                                       ind_imag_skin_metrics, one_way_summary, local_cli.coil_source_scale,
                                       local_cli.level == Team7ValidationLevel::CoilOnly
                                           ? 0.0
                                           : hostOphelieEdgeFluxFallbackFraction(plate_body.getBaseParticles(),
                                                                                 plate_names));
        if (!x_band_metrics.empty())
        {
            appendTeam7ValidationXBandRecord(history_path, x_band_metrics);
        }
        writeTeam7SummaryText(team7OutputArtifactPath("team7_summary_", local_cli.level, output_tag, ".txt"), level_tag,
                              output_tag, pass_report, coil_metrics, total_phase90_metrics, one_way_summary,
                              local_cli.coil_source_scale);
        appendTeam7ValidationRunRecord("./output/team7_validation_history.txt", level_tag, output_tag, pass_report,
                                       coil_metrics, coil_span_metrics, total_metrics, total_phase90_metrics,
                                       total_magnitude_metrics, ind_imag_skin_metrics, one_way_summary,
                                       local_cli.coil_source_scale,
                                       local_cli.level == Team7ValidationLevel::CoilOnly
                                           ? 0.0
                                           : hostOphelieEdgeFluxFallbackFraction(plate_body.getBaseParticles(),
                                                                                 plate_names));
        if (!x_band_metrics.empty())
        {
            appendTeam7ValidationXBandRecord("./output/team7_validation_history.txt", x_band_metrics);
        }
    }

    if (local_cli.level != Team7ValidationLevel::CoilOnly)
    {
        BodyStatesRecordingToVtp plate_vtp(sph_system);
        plate_vtp.addToWrite<Vecd>(plate_body, plate_names.b_coil_real);
        plate_vtp.addToWrite<Vecd>(plate_body, plate_names.b_coil_imag);
        plate_vtp.addToWrite<Vecd>(plate_body, plate_names.b_ind_real);
        plate_vtp.addToWrite<Vecd>(plate_body, plate_names.b_ind_imag);
        plate_vtp.addToWrite<Vecd>(plate_body, plate_names.j_imag);
        plate_vtp.addToWrite<Vecd>(plate_body, plate_names.j_edge_recon_imag);
        plate_vtp.addToWrite<Real>(plate_body, plate_names.phi_imag);
        plate_vtp.addToWrite<Real>(plate_body, kTeam7DepthBinZMm);
        plate_body.setNewlyUpdated();
        plate_vtp.writeToFile(0);
        logOphelieOutputArtifact(IO::getEnvironment().OutputFolder() + "/PlateBody_ite_0000000000.vtp");
    }
    return pass_report.smoke_passed ? 0 : 1;
}
