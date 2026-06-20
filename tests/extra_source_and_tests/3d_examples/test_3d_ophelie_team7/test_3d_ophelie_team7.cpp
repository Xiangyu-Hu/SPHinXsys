/**
 * @file test_3d_ophelie_team7.cpp
 * @brief OPHELIE-like induction on TEAM7-like geometry: conducting plate + annular coil + air domain.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_team7_geometry.h"
#include "electromagnetic_ophelie_team7_coil_path_source.h"
#include "electromagnetic_ophelie_team7_native_geometry.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"
#include "electromagnetic_ophelie_team7_probe.h"
#include "electromagnetic_ophelie_racetrack_source.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_relaxation.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <iostream>
#include <memory>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

inline OphelieRunMetrics collectPlateMetrics(BaseParticles &plate_particles, const OphelieGlassFieldNames &names)
{
    const size_t n = plate_particles.TotalRealParticles();
    syncGlassElectromagneticFieldsToHost(plate_particles, names);

    OphelieRunMetrics metrics;
    metrics.n_glass = n;
    metrics.max_a_src = hostVecdFieldMax(plate_particles, names.a_src_real, n);
    metrics.max_b_src = hostVecdFieldMax(plate_particles, names.b_src_real, n);
    metrics.max_e_imag = hostVecdFieldMax(plate_particles, names.e_imag, n);
    metrics.max_j_imag = hostVecdFieldMax(plate_particles, names.j_imag, n);
    metrics.max_joule_heat = hostScalarFieldMax(plate_particles, names.joule_heat, n);

    syncVariableToHost<Real>(plate_particles, names.joule_heat);
    const Real *joule = plate_particles.getVariableDataByName<Real>(names.joule_heat);
    metrics.min_joule_heat = joule[0];
    for (size_t i = 1; i != n; ++i)
    {
        metrics.min_joule_heat = std::min(metrics.min_joule_heat, joule[i]);
    }
    return metrics;
}

inline void applyTeam7DefaultParameters(OphelieParameters &params, const OphelieTeam7Geometry &geom)
{
    params.coil_center_ = geom.coil_center_;
    params.glass_center_ = geom.plate_center_;
    params.frequency_ = 50.0;
    /** Physical Al ~3.54e7 S/m; SPH phi solve needs coarser sigma or smaller dp — override with --sigma=. */
    params.sigma_glass_ = 1.0e4;
    params.current_amplitude_ = 1.0;
    params.number_of_turns_ = 1.0;
    params.coil_j0_override_ =
        params.number_of_turns_ * params.current_amplitude_ / (geom.coilCurrentCrossSectionArea() + TinyReal);
    params.target_joule_power_ = 50.0e3;
    params.enable_phi_correction_ = true;
}

} // namespace

int main(int ac, char *av[])
{
    logOphelieRunContext();

    OphelieTeam7Geometry analytic_geom;
    OphelieParameters params;
    applyTeam7DefaultParameters(params, analytic_geom);

    OphelieTestCliOptions cli_options;
    const StdVec<std::string> filtered_arguments = filterOphelieTestCommandLine(ac, av, params, cli_options);
    if (cli_options.no_power_scaling)
    {
        params.enable_power_scaling_ = false;
    }
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
    native_mesh.use_small_air_box_ = cli_options.native_small_air_box;
    native_mesh.relax_air_particles_ = cli_options.native_relax_air;

    const size_t relaxation_steps =
        cli_options.relaxation_steps > 0 ? cli_options.relaxation_steps
                                         : (cli_options.native_stl ? native_mesh.relaxation_steps_ : 400);
    const size_t relaxation_log_every =
        cli_options.relaxation_log_every > 0 ? cli_options.relaxation_log_every : 50;
    const size_t relaxation_vtp_every = cli_options.relaxation_vtp_every;
    native_mesh.relaxation_steps_ = relaxation_steps;
    native_mesh.relaxation_log_every_ = relaxation_log_every;
    native_mesh.relaxation_vtp_every_ = relaxation_vtp_every;
    if (cli_options.team7_coil_turns > TinyReal)
    {
        native_mesh.team7_coil_turns_ = cli_options.team7_coil_turns;
    }

    Real dp = 0.04;
    BoundingBoxd system_bounds(Vecd::Zero(), Vecd::Ones());
    if (cli_options.native_stl)
    {
        const Vec3d air_lower =
            native_mesh.use_small_air_box_ ? native_mesh.air_lower_small_ : native_mesh.air_lower_standard_;
        const Vec3d air_upper =
            native_mesh.use_small_air_box_ ? native_mesh.air_upper_small_ : native_mesh.air_upper_standard_;
        system_bounds = BoundingBoxd(native_mesh.stl_scale_to_meter_ * air_lower,
                                       native_mesh.stl_scale_to_meter_ * air_upper);
        dp = native_mesh.dp_ref_mm_ * native_mesh.stl_scale_to_meter_;
        std::cout << "[ophelie] TEAM7 native STL mode dp=" << dp << " m (dp_mm=" << native_mesh.dp_ref_mm_
                  << ") coil_stl=" << team7NativeCoilStlPath() << " plate_stl=" << team7NativePlateStlPath()
                  << std::endl;
    }
    else
    {
        const Real boundary_width = 3.0 * dp;
        system_bounds = BoundingBoxd(analytic_geom.domainMin(boundary_width), analytic_geom.domainMax(boundary_width));
    }

    SPHSystem sph_system(system_bounds, dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    if (!cli_options.reload_dir.empty())
    {
        IO::getEnvironment().resetReloadFolder(cli_options.reload_dir, true);
        std::cout << "[ophelie] reload folder: " << IO::getEnvironment().ReloadFolder() << std::endl;
    }
    sph_system.handleCommandlineOptions(filtered_ac, filtered_av);

    const std::string coil_stl_path = team7NativeCoilStlPath();
    const std::string plate_stl_path = team7NativePlateStlPath();
    const Real coil_adapt_ratio =
        cli_options.native_stl ? (native_mesh.dp_ref_mm_ / native_mesh.dp_coil_mm_) : 1.0;
    const Real plate_adapt_ratio =
        cli_options.native_stl ? (native_mesh.dp_ref_mm_ / native_mesh.dp_plate_mm_) : 1.0;

    OphelieProgressLogger setup_progress("setup");
    SharedPtr<ComplexShape> plate_shape;
    SharedPtr<ComplexShape> coil_shape;
    if (cli_options.native_stl)
    {
        plate_shape = makeShared<OphelieTeam7NativePlateShape>("PlateBody", plate_stl_path, native_mesh.stl_scale_to_meter_);
        coil_shape = makeShared<OphelieTeam7NativeCoilShape>("CoilSourceBody", coil_stl_path, native_mesh.stl_scale_to_meter_);
    }
    else
    {
        plate_shape = makeShared<OphelieTeam7PlateShape>("PlateBody", analytic_geom);
        coil_shape = makeShared<OphelieTeam7AnnularCoilShape>("CoilSourceBody", analytic_geom);
    }

    SolidBody plate_body(sph_system, plate_shape);
    plate_body.defineAdaptation<SPHAdaptation>(1.15, plate_adapt_ratio);
    plate_body.defineMatterMaterial<Solid>();
#if SPHINXSYS_USE_SYCL
    LevelSetShape &plate_level_set =
        plate_body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
#else
    LevelSetShape &plate_level_set = plate_body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
#endif

    SolidBody coil_body(sph_system, coil_shape);
    coil_body.defineAdaptation<SPHAdaptation>(1.15, coil_adapt_ratio);
    coil_body.defineMatterMaterial<Solid>();
#if SPHINXSYS_USE_SYCL
    LevelSetShape &coil_level_set =
        coil_body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
#else
    LevelSetShape &coil_level_set = coil_body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
#endif

    std::unique_ptr<RealBody> air_body;
    LevelSetShape *air_level_set = nullptr;
    /** Air particles are only for native SYCL relax; do not register AirBody on --reload EM runs (avoids segfault). */
    if (cli_options.native_stl && native_mesh.relax_air_particles_ && sph_system.RunParticleRelaxation())
    {
        const Vec3d air_lower =
            native_mesh.use_small_air_box_ ? native_mesh.air_lower_small_ : native_mesh.air_lower_standard_;
        const Vec3d air_upper =
            native_mesh.use_small_air_box_ ? native_mesh.air_upper_small_ : native_mesh.air_upper_standard_;
        air_body = std::make_unique<RealBody>(
            sph_system, makeShared<OphelieTeam7NativeAirShape>("AirBody", air_lower, air_upper, coil_stl_path,
                                                               plate_stl_path, native_mesh.stl_scale_to_meter_));
#if SPHINXSYS_USE_SYCL
        air_level_set = &air_body->defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
#else
        air_level_set = &air_body->defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
#endif
    }

    Vecd coil_center = analytic_geom.coil_center_;
    Real div_j_characteristic_length = analytic_geom.plate_radius_;
    OphelieTeam7NativeDerivedGeometry native_derived;
    OphelieTeam7CoilPathPrepareSummary coil_path_summary;
    OphelieTeam7CoilPathSourceSpec coil_path_spec;
    bool coil_path_audit_ok = true;
    const bool use_particle_reload = !sph_system.RunParticleRelaxation() && sph_system.ReloadParticles();

    if (sph_system.RunParticleRelaxation())
    {
        plate_body.generateParticles<BaseParticles, Lattice>();
        coil_body.generateParticles<BaseParticles, Lattice>();
        if (air_body)
        {
            air_body->generateParticles<BaseParticles, Lattice>();
        }
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();

        if (cli_options.native_stl)
        {
#if SPHINXSYS_USE_SYCL
            relaxTeam7NativeStlBodies(sph_system, coil_body, plate_body, coil_level_set, plate_level_set, native_mesh,
                                      air_body.get(), air_level_set);
#else
            relaxSolidBodyParticles(sph_system, plate_body, plate_level_set, "PlateBody", relaxation_steps,
                                    relaxation_log_every, relaxation_vtp_every);
            relaxSolidBodyParticles(sph_system, coil_body, coil_level_set, "CoilSourceBody", relaxation_steps,
                                    relaxation_log_every, relaxation_vtp_every);
#endif
        }
        else
        {
            relaxSolidBodyParticles(sph_system, plate_body, plate_level_set, "PlateBody", relaxation_steps,
                                    relaxation_log_every, relaxation_vtp_every);
            relaxSolidBodyParticles(sph_system, coil_body, coil_level_set, "CoilSourceBody", relaxation_steps,
                                    relaxation_log_every, relaxation_vtp_every);
        }
        computeTeam7CoilPlateNormalsFromShape(sph_system, coil_body, plate_body);
        ReloadParticleIO write_reload({&coil_body, &plate_body});
        registerTeam7CoilPlateNormalsForReload(write_reload, coil_body, plate_body);
        write_reload.writeToFile(0);
        std::cout << "test_3d_ophelie_team7 particle relaxation finished. Reload.xml (CoilSourceBody+PlateBody+NormalDirection) -> "
                  << IO::getEnvironment().ReloadFolder() << std::endl;
        return 0;
    }

    if (use_particle_reload)
    {
        reloadTeam7CoilPlateParticlesWithNormals(coil_body, plate_body);
        std::cout << "[ophelie] loaded relaxed particles from " << IO::getEnvironment().ReloadFolder()
                  << "/Reload.xml" << std::endl;
    }
    else
    {
        plate_body.generateParticles<BaseParticles, Lattice>();
        coil_body.generateParticles<BaseParticles, Lattice>();
        if (air_body)
        {
            air_body->generateParticles<BaseParticles, Lattice>();
        }
        if (!cli_options.skip_relaxation)
        {
            sph_system.initializeSystemCellLinkedLists();
            sph_system.initializeSystemConfigurations();
            if (cli_options.native_stl)
            {
#if SPHINXSYS_USE_SYCL
                relaxTeam7NativeStlBodies(sph_system, coil_body, plate_body, coil_level_set, plate_level_set,
                                          native_mesh, air_body.get(), air_level_set);
#else
                relaxSolidBodyParticles(sph_system, plate_body, plate_level_set, "PlateBody", relaxation_steps,
                                        relaxation_log_every, relaxation_vtp_every);
                relaxSolidBodyParticles(sph_system, coil_body, coil_level_set, "CoilSourceBody", relaxation_steps,
                                        relaxation_log_every, relaxation_vtp_every);
#endif
            }
            else
            {
                relaxSolidBodyParticles(sph_system, plate_body, plate_level_set, "PlateBody", relaxation_steps,
                                        relaxation_log_every, relaxation_vtp_every);
                relaxSolidBodyParticles(sph_system, coil_body, coil_level_set, "CoilSourceBody", relaxation_steps,
                                        relaxation_log_every, relaxation_vtp_every);
            }
        }
        else
        {
            std::cout << "[ophelie] --skip-relax: using lattice particles without relaxation" << std::endl;
        }
    }

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    if (cli_options.native_stl)
    {
        native_derived = deriveTeam7NativeGeometry(coil_body, plate_body);
        applyTeam7NativeParameters(params, native_mesh, native_derived, coil_body.getBaseParticles());
        if (cli_options.coil_source_model == OphelieCoilSourceModel::VolumeRacetrack)
        {
            coil_path_spec.turns = native_mesh.team7_coil_turns_;
            coil_path_spec.current_per_turn = native_mesh.team7_coil_current_per_turn_;
            coil_path_summary =
                prepareOphelieTeam7VolumeRacetrackCoilSource(coil_body, native_derived.coil_bbox_, coil_path_spec);
            params.coil_j0_override_ = coil_path_summary.j0;
            printOphelieTeam7CoilPathPrepareSummary(coil_path_summary, coil_path_spec);
        }
        if (!cli_options.sigma_user_set)
        {
            params.sigma_glass_ = 3.54e7;
        }
        std::cout << "[ophelie] native TEAM7 coil_source_model="
                  << ophelieCoilSourceModelName(cli_options.coil_source_model)
                  << " coil_center=" << native_derived.coil_center_.transpose()
                  << " A_cross=" << native_derived.coil_current_cross_section_m2_
                  << " coil_volume=" << native_derived.coil_volume_m3_
                  << " mean_radius=" << native_derived.coil_mean_radius_m_
                  << " max_xy_radius=" << native_derived.coil_max_xy_radius_m_
                  << " J0=" << params.coil_j0_override_ << " sigma=" << params.sigma_glass_ << std::endl;
        coil_center = native_derived.coil_center_;
        div_j_characteristic_length = native_derived.plate_characteristic_length_m_;
        if (cli_options.compare_team7_bz)
        {
            printTeam7NativeStlMeshBBoxAudit(Team7NativeStlMeshBBoxMm{});
            printTeam7NativeGeometryProbeAudit(coil_body, plate_body, native_derived, params, native_mesh);
        }
    }

    setup_progress.log(std::string(cli_options.native_stl ? "native-stl " : "analytic ") + "n_plate=" +
                       std::to_string(plate_body.getBaseParticles().TotalRealParticles()) + " n_coil=" +
                       std::to_string(coil_body.getBaseParticles().TotalRealParticles()) + " dp=" + std::to_string(dp) +
                       (use_particle_reload ? " source=reload" : ""));
    setup_progress.finish();

    OphelieCoilFieldNames coil_names;
    OphelieGlassFieldNames plate_names;

    RegisterOphelieCoilFields register_coil_fields(coil_body, coil_names);
    RegisterOphelieGlassFields register_plate_fields(plate_body, plate_names);
    (void)register_coil_fields;
    (void)register_plate_fields;

    UniquePtr<Inner<>> plate_inner = makeUnique<Inner<>>(plate_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_plate_sigma(plate_body, plate_names,
                                                                                     params.sigma_glass_);
    StateDynamics<MainExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot_savart(
        plate_body, coil_body, plate_names, coil_names, params);

    syncCoilSourceFieldsToDevice(coil_body.getBaseParticles(), coil_names);
    syncGlassElectromagneticFieldsToDevice(plate_body.getBaseParticles(), plate_names);

    Team7BzCompareMetrics team7_bz_metrics;
    team7_bz_metrics.passed = true;

    {
        OphelieProgressLogger em_progress("em_init");
        assign_plate_sigma.exec();
        if (cli_options.coil_source_model == OphelieCoilSourceModel::VolumeRacetrack)
        {
            StateDynamics<MainExecutionPolicy, InitializeOphelieVolumeRacetrackCoilSourceCK> initialize_volume_racetrack(
                coil_body, coil_names, coil_path_summary.j0);
            initialize_volume_racetrack.exec();
            const Real integrated_current =
                hostCoilIntegratedCurrentFromJSrc(coil_body.getBaseParticles(), coil_names, coil_path_summary.path_length_m);
            coil_path_audit_ok =
                ophelieTeam7CoilPathAmpereTurnsAuditPassed(coil_path_summary, integrated_current);
        }
        else
        {
            StateDynamics<MainExecutionPolicy, InitializeOphelieCoilSourceCK> initialize_coil_source(
                coil_body, coil_names, params, coil_center);
            initialize_coil_source.exec();
        }
        em_progress.log("coil/plate source fields assigned");
        compute_biot_savart.exec();
        em_progress.log("biot savart coil->plate done");
        em_progress.finish();
    }

    if (cli_options.compare_team7_bz)
    {
        if (!cli_options.native_stl)
        {
            std::cout << "[ophelie] WARNING: --compare-team7-bz is intended for --native-stl geometry." << std::endl;
        }
        if (params.enable_phi_correction_)
        {
            std::cout << "[ophelie] WARNING: Bz probe uses coil-only Biot; use --no-phi for reference compare." << std::endl;
        }
        const std::string reference_dir =
            cli_options.team7_reference_dir.empty()
                ? "tests/extra_source_and_tests/3d_examples/reference_data/team7"
                : cli_options.team7_reference_dir;
        StdVec<Team7BzProbePoint> reference_probes;
        if (loadTeam7BzA1B1Reference(reference_dir, Team7ReferenceProbeLineMm::y_mm, Team7ReferenceProbeLineMm::z_mm,
                                     1.0e-3, reference_probes))
        {
            const Real total_current =
                native_mesh.team7_coil_turns_ * native_mesh.team7_coil_current_per_turn_;
            const Team7NativeStlMeshBBoxMm mesh_bbox;

            if (cli_options.coil_source_model == OphelieCoilSourceModel::FilamentRacetrack)
            {
                if (cli_options.racetrack_sweep || cli_options.compare_team7_bz_rect_loop)
                {
                    runFilamentRacetrackTeam7BzSweep(mesh_bbox, total_current, params.mu0_, params.softening_length_,
                                                       reference_probes, cli_options.racetrack_ds_mm,
                                                       cli_options.team7_bz_rms_smoke_threshold);
                    team7_bz_metrics.passed = true;
                }
                else
                {
                    OphelieRacetrackParams racetrack;
                    racetrack.inset_mm = cli_options.racetrack_inset_mm;
                    racetrack.z_mm = cli_options.racetrack_z_mm;
                    racetrack.ds_mm = cli_options.racetrack_ds_mm;
                    const Team7BzDiagnosticSummary summary = evaluateFilamentRacetrackTeam7BzDiagnostic(
                        mesh_bbox, racetrack, total_current, params.mu0_, params.softening_length_, reference_probes,
                        cli_options.team7_bz_rms_smoke_threshold);
                    printTeam7BzDiagnosticSummary(summary);
                    team7_bz_metrics = summary.metrics;
                }
            }
            else
            {
                StdVec<Team7BzProbePoint> simulated_probes;
                evaluateCoilBiotSavartBzAtProbes(coil_body.getBaseParticles(), coil_names, params, reference_probes,
                                                  simulated_probes);
                team7_bz_metrics = compareTeam7BzPhase0(simulated_probes, cli_options.team7_bz_rms_smoke_threshold);
                OphelieRacetrackParams racetrack;
                const Team7BzDiagnosticSummary summary =
                    summarizeTeam7BzDiagnostic(simulated_probes, team7_bz_metrics, cli_options.coil_source_model,
                                               racetrack, cli_options.team7_bz_rms_smoke_threshold);
                printTeam7BzDiagnosticSummary(summary);
                printTeam7BzCompareReport(simulated_probes, team7_bz_metrics, cli_options.team7_bz_rms_smoke_threshold);
                if (cli_options.native_stl)
                {
                    const Real coil_x_min_mm = native_derived.coil_bbox_.lower_[0] * 1000.0;
                    const Team7BzCompareMetrics coil_span_metrics = compareTeam7BzPhase0OverXRange(
                        simulated_probes, coil_x_min_mm, Team7ReferenceProbeLineMm::x_end_mm,
                        cli_options.team7_bz_rms_smoke_threshold);
                    std::cout << "[ophelie] TEAM7 Bz subset x_mm>=" << coil_x_min_mm
                              << " (under coil x-span) n=" << coil_span_metrics.n_probes
                              << " rms_rel_err=" << coil_span_metrics.rms_rel_error
                              << " passed=" << (coil_span_metrics.passed ? 1 : 0) << std::endl;
                }
                writeTeam7BzCompareCsv("./output/team7_bz_A1B1_compare.csv", simulated_probes);
                StdVec<Team7PlateBzNearProbeSample> plate_samples;
                sampleNearestPlateCoilBzAtProbes(plate_body.getBaseParticles(), plate_names, simulated_probes, 1.0e-3,
                                                 plate_samples);
                printTeam7PlateBzNearProbeReport(plate_samples);
            }

            if (cli_options.compare_team7_bz_loop && cli_options.native_stl &&
                cli_options.coil_source_model == OphelieCoilSourceModel::VolumeETheta)
            {
                StdVec<Team7BzProbePoint> loop_probes;
                evaluateTeam7CircularLoopBiotSavartBzAtProbes(native_derived.coil_center_,
                                                                native_derived.coil_mean_radius_m_, total_current,
                                                                params.mu0_, reference_probes, loop_probes);
                const Team7BzCompareMetrics loop_metrics =
                    compareTeam7BzPhase0(loop_probes, cli_options.team7_bz_rms_smoke_threshold);
                std::cout << "[ophelie] TEAM7 Bz circular-loop diagnostic (R=mean_radius, I=N*I_turn) n="
                          << loop_metrics.n_probes << " rms_rel_err=" << loop_metrics.rms_rel_error
                          << " peak_ref_x=" << reference_probes[loop_metrics.peak_ref_index].x_mm
                          << " peak_loop_x=" << loop_probes[loop_metrics.peak_sim_index].x_mm << std::endl;
                writeTeam7BzCompareCsv("./output/team7_bz_A1B1_loop_compare.csv", loop_probes);
            }
        }
        else
        {
            team7_bz_metrics.passed = false;
        }
    }

    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        plate_body, plate_names);
    combine_vector_potential.exec();

    OphelieRunMetrics metrics;
    Real phi_solver_rel_residual = 0.0;
    Real self_induction_j_rel_change = 0.0;
    size_t self_induction_iterations_used = 0;

    BaseParticles &plate_particles = plate_body.getBaseParticles();
    const size_t n_plate = plate_particles.TotalRealParticles();
    const size_t n_coil = coil_body.getBaseParticles().TotalRealParticles();

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> compute_ejq_no_phi(plate_body, plate_names, params);
    OphelieDivJMetrics div_j_level0_metrics;
    OphelieDivJMetrics div_j_phi_metrics;

    if (params.enable_self_induction_)
    {
        Real self_induction_phi_eq_res_vol = 0.0;
        bool self_induction_picard_converged = false;
        self_induction_j_rel_change = runOphelieSelfInductionWithPhiSolve<MainExecutionPolicy>(
            plate_body, *plate_inner, plate_names, params, phi_solver_rel_residual, self_induction_iterations_used,
            self_induction_phi_eq_res_vol, self_induction_picard_converged);
        div_j_phi_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(plate_body, *plate_inner, plate_names, div_j_characteristic_length);
    }
    else if (params.enable_phi_correction_)
    {
        compute_ejq_no_phi.exec();
        div_j_level0_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(plate_body, *plate_inner, plate_names, div_j_characteristic_length);

        OphelieProgressLogger phi_progress("phi_solve");
        phi_solver_rel_residual = solvePhiImag<MainExecutionPolicy>(plate_body, *plate_inner, plate_names, params);
        phi_progress.finish("rel_res=" + std::to_string(phi_solver_rel_residual));

        InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
            *plate_inner, plate_names);
        StateDynamics<MainExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(plate_body, plate_names, params);

        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(plate_body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*plate_inner);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        compute_grad_phi.exec();
        compute_ejq_with_phi.exec();
        div_j_phi_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(plate_body, *plate_inner, plate_names, div_j_characteristic_length);
    }
    else
    {
        compute_ejq_no_phi.exec();
        div_j_level0_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(plate_body, *plate_inner, plate_names, div_j_characteristic_length);
    }

    const Real joule_power_raw = hostVolWeightedSum(plate_particles, plate_names.joule_heat, n_plate);
    const OpheliePowerScalingFactors scaling_factors = computeOpheliePowerScalingFactors(params, joule_power_raw);
    if (params.enable_power_scaling_)
    {
        StateDynamics<MainExecutionPolicy, ScaleOphelieElectromagneticFieldsCK> scale_em_fields(
            plate_body, plate_names, scaling_factors.field_scale, scaling_factors.power_scale);
        scale_em_fields.exec();
    }
    else
    {
        std::cout << "[ophelie] power_scaling disabled (raw fields, I_eff=I0)" << std::endl;
    }

    metrics = collectPlateMetrics(plate_particles, plate_names);
    metrics.n_coil = n_coil;
    metrics.joule_power_raw = joule_power_raw;
    metrics.joule_power_scaled = hostVolWeightedSum(plate_particles, plate_names.joule_heat, n_plate);
    metrics.power_scale = scaling_factors.power_scale;
    metrics.field_scale = scaling_factors.field_scale;
    metrics.effective_current_amplitude = scaling_factors.effective_current_amplitude;
    metrics.div_j_rel_level0 = div_j_level0_metrics.div_j_rel;
    metrics.div_j_rel_phi = div_j_phi_metrics.div_j_rel > 0.0 ? div_j_phi_metrics.div_j_rel : div_j_level0_metrics.div_j_rel;
    if (metrics.div_j_rel_phi > TinyReal && metrics.div_j_rel_level0 > TinyReal)
    {
        metrics.div_j_reduction = metrics.div_j_rel_level0 / metrics.div_j_rel_phi;
    }
    std::cout << "[ophelie] power_scaling P_raw=" << metrics.joule_power_raw << " P_target=" << params.target_joule_power_
              << " power_scale=" << metrics.power_scale << " field_scale=" << metrics.field_scale
              << " I_eff=" << metrics.effective_current_amplitude << std::endl;
    std::cout << "[ophelie] divJ_rel_level0=" << metrics.div_j_rel_level0 << " divJ_rel_phi=" << metrics.div_j_rel_phi
              << " divJ_reduction=" << metrics.div_j_reduction << std::endl;
    metrics.with_phi_correction = params.enable_phi_correction_;
    metrics.max_a_coil = hostVecdFieldMax(plate_particles, plate_names.a_coil_real, n_plate);
    metrics.max_a_ind = hostVecdFieldMax(plate_particles, plate_names.a_ind_real, n_plate);
    if (params.enable_phi_correction_ || params.enable_self_induction_)
    {
        metrics.phi_solver_rel_residual = phi_solver_rel_residual;
        metrics.phi_solver_kind = params.phi_solver_kind_;
        metrics.self_induction_j_rel_change = self_induction_j_rel_change;
        metrics.self_induction_iterations_used = self_induction_iterations_used;
        metrics.max_phi_imag = hostScalarFieldMax(plate_particles, plate_names.phi_imag, n_plate);
        metrics.max_phi_rhs_imag = hostScalarFieldMax(plate_particles, plate_names.phi_rhs_imag, n_plate);
    }

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(coil_body, coil_names.j_src_real);
    write_states.addToWrite<Vecd>(plate_body, plate_names.a_src_real);
    write_states.addToWrite<Vecd>(plate_body, plate_names.b_src_real);
    write_states.addToWrite<Vecd>(plate_body, plate_names.e_imag);
    write_states.addToWrite<Vecd>(plate_body, plate_names.j_imag);
    write_states.addToWrite<Real>(plate_body, plate_names.joule_heat);
    write_states.addToWrite<Real>(plate_body, plate_names.sigma);
    if (params.enable_phi_correction_)
    {
        write_states.addToWrite<Real>(plate_body, plate_names.phi_imag);
        write_states.addToWrite<Vecd>(plate_body, plate_names.grad_phi_imag);
    }
    {
        OphelieProgressLogger vtp_progress("vtp_output");
        write_states.writeToFile(0);
        vtp_progress.finish("wrote ./output (see cwd)");
    }

    const bool phi_residual_ok =
        !params.enable_phi_correction_ || cli_options.ophelie_smoke ||
        (metrics.max_phi_imag > 0.0 && metrics.max_phi_rhs_imag > 0.0 &&
         metrics.phi_solver_rel_residual <
             10.0 * (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES
                         ? params.phi_gmres_tolerance_
                         : params.phi_solver_kind_ == OpheliePhiSolverKind::PCG
                               ? params.phi_pcg_tolerance_
                               : params.phi_jacobi_tolerance_) &&
         (!params.enable_self_induction_ ||
          metrics.self_induction_j_rel_change < params.self_induction_j_tolerance_));

    const bool team7_bz_ok = !cli_options.compare_team7_bz || team7_bz_metrics.passed;
    const bool coil_path_ok = cli_options.coil_source_model != OphelieCoilSourceModel::VolumeRacetrack ||
                              !cli_options.native_stl || coil_path_audit_ok;
    const bool power_ok = !params.enable_power_scaling_ || (std::isfinite(metrics.joule_power_scaled) &&
                                                           metrics.joule_power_scaled > 0.0);
    const bool passed = metrics.n_glass > 0 && metrics.n_coil > 0 && metrics.max_a_src > 0.0 && metrics.max_b_src > 0.0 &&
                        metrics.max_e_imag > 0.0 && metrics.max_j_imag > 0.0 && metrics.max_joule_heat > 0.0 &&
                        metrics.min_joule_heat >= 0.0 && power_ok && phi_residual_ok && team7_bz_ok && coil_path_ok;

    std::cout << "test_3d_ophelie_team7"
              << (cli_options.native_stl ? " native_stl=1" : " native_stl=0")
              << " coil_source_model=" << ophelieCoilSourceModelName(cli_options.coil_source_model) << " dp=" << dp
              << " n_plate=" << metrics.n_glass << " n_coil=" << metrics.n_coil
              << " frequency=" << params.frequency_ << " sigma=" << params.sigma_glass_
              << " J0=" << params.coil_j0_override_ << " phi_correction=" << (params.enable_phi_correction_ ? 1 : 0)
              << " phi_solver=" << phiSolverKindName(params.phi_solver_kind_)
              << " phi_rel_res=" << metrics.phi_solver_rel_residual
              << " self_ind_iters=" << metrics.self_induction_iterations_used
              << " self_ind_J_rel=" << metrics.self_induction_j_rel_change << " max_ACoil=" << metrics.max_a_coil
              << " max_AInd=" << metrics.max_a_ind << " max_PhiImag=" << metrics.max_phi_imag
              << " max_PhiRhs=" << metrics.max_phi_rhs_imag << " max_ASrc=" << metrics.max_a_src
              << " max_BSrc=" << metrics.max_b_src << " max_EImag=" << metrics.max_e_imag << " max_JImag=" << metrics.max_j_imag
              << " P_raw=" << metrics.joule_power_raw << " P_scaled=" << metrics.joule_power_scaled
              << " power_scale=" << metrics.power_scale << " field_scale=" << metrics.field_scale
              << " I_eff=" << metrics.effective_current_amplitude << " divJ_L0=" << metrics.div_j_rel_level0
              << " divJ_phi=" << metrics.div_j_rel_phi << " divJ_red=" << metrics.div_j_reduction
              << " target_P=" << params.target_joule_power_ << " min_Joule=" << metrics.min_joule_heat
              << " max_Joule=" << metrics.max_joule_heat
              << " team7_bz_rms=" << team7_bz_metrics.rms_rel_error
              << " team7_bz_passed=" << (team7_bz_metrics.passed ? 1 : 0)
              << " coil_path_audit=" << (coil_path_audit_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;

    return passed ? 0 : 1;
}
