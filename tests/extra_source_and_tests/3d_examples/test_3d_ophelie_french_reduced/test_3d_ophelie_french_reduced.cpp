/**
 * @file test_3d_ophelie_french_reduced.cpp
 * @brief French-paper-inspired reduced case: mesh glass cylinder + multiloop line-source coil.
 *
 * Reference: Jacoutot et al., Chem. Eng. Process. 47 (2008) 449-455.
 * This is NOT an exact CAD reconstruction; see docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_glass_mesh_relax.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_progress.h"
#include "electromagnetic_ophelie_relaxation.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

inline OphelieRunMetrics collectGlassMetrics(BaseParticles &glass_particles, const OphelieGlassFieldNames &glass_names)
{
    const size_t n_glass = glass_particles.TotalRealParticles();
    syncGlassElectromagneticFieldsToHost(glass_particles, glass_names);

    OphelieRunMetrics metrics;
    metrics.n_glass = n_glass;
    metrics.max_a_src = hostVecdFieldMax(glass_particles, glass_names.a_src_real, n_glass);
    metrics.max_b_src = hostVecdFieldMax(glass_particles, glass_names.b_src_real, n_glass);
    metrics.max_e_imag = hostVecdFieldMax(glass_particles, glass_names.e_imag, n_glass);
    metrics.max_j_imag = hostVecdFieldMax(glass_particles, glass_names.j_imag, n_glass);
    metrics.max_joule_heat = hostScalarFieldMax(glass_particles, glass_names.joule_heat, n_glass);

    syncVariableToHost<Real>(glass_particles, glass_names.joule_heat);
    const Real *joule = glass_particles.getVariableDataByName<Real>(glass_names.joule_heat);
    metrics.min_joule_heat = joule[0];
    for (size_t i = 1; i != n_glass; ++i)
    {
        metrics.min_joule_heat = std::min(metrics.min_joule_heat, joule[i]);
    }
    return metrics;
}

inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return defineFrenchGlassMeshLevelSet(body, false);
}

inline void relaxFrenchGlassMeshBody(SPHSystem &sph_system, SolidBody &glass_body, LevelSetShape &glass_level_set,
                                     size_t relaxation_steps, size_t relaxation_vtp_every)
{
#if SPHINXSYS_USE_SYCL
    runFrenchGlassSphinxsysStyleRelax(sph_system, glass_body, glass_level_set, relaxation_steps, relaxation_vtp_every);
#else
    (void)sph_system;
    (void)glass_body;
    (void)glass_level_set;
    (void)relaxation_steps;
    (void)relaxation_vtp_every;
    std::cout << "[ophelie] ERROR: mesh-cylinder relax requires SYCL." << std::endl;
    std::exit(1);
#endif
}

inline void generateFrenchReducedVisualLattice(SolidBody *coil_visual_body, SolidBody *crucible_visual_body)
{
    if (coil_visual_body)
    {
        coil_visual_body->generateParticles<BaseParticles, Lattice>();
    }
    if (crucible_visual_body)
    {
        crucible_visual_body->generateParticles<BaseParticles, Lattice>();
    }
}

inline void relaxFrenchReducedBodies(SPHSystem &sph_system, SolidBody &glass_body, LevelSetShape &glass_level_set,
                                     SolidBody *coil_visual_body, LevelSetShape *coil_visual_level_set,
                                     SolidBody *crucible_visual_body, LevelSetShape *crucible_visual_level_set,
                                     size_t relaxation_steps, size_t relaxation_log_every, size_t relaxation_vtp_every)
{
    relaxFrenchGlassMeshBody(sph_system, glass_body, glass_level_set, relaxation_steps, relaxation_vtp_every);
    if (coil_visual_body && coil_visual_level_set)
    {
        relaxSolidBodyParticles(sph_system, *coil_visual_body, *coil_visual_level_set, "CoilVisualBody",
                                relaxation_steps, relaxation_log_every, relaxation_vtp_every);
    }
    if (crucible_visual_body && crucible_visual_level_set)
    {
        relaxSolidBodyParticles(sph_system, *crucible_visual_body, *crucible_visual_level_set, "CrucibleWallVisualBody",
                                relaxation_steps, relaxation_log_every, relaxation_vtp_every);
    }
}

inline void writeFrenchReducedReload(SolidBody &glass_body, SolidBody *coil_visual_body,
                                     SolidBody *crucible_visual_body)
{
    SPHBodyVector bodies{&glass_body};
    if (coil_visual_body)
    {
        bodies.push_back(coil_visual_body);
    }
    if (crucible_visual_body)
    {
        bodies.push_back(crucible_visual_body);
    }
    ReloadParticleIO write_reload(bodies);
    write_reload.writeToFile(0);
}

} // namespace

int main(int ac, char *av[])
{
    logOphelieRunContext();

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);

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
    const StdVec<std::string> filtered_arguments =
        filterOphelieTestCommandLine(french_ac, french_av.data(), params, cli_options);
    if (cli_options.no_power_scaling)
    {
        params.enable_power_scaling_ = false;
    }

    OphelieFrenchLiteratureProfile literature_profile;
    literature_profile.calibrate_coil_current = cli_options.literature_calibrate_current;
    literature_profile.div_j_l2_reduction_min = cli_options.literature_div_j_l2_reduction_min;
    if (cli_options.literature_mode)
    {
        applyFrenchLiteratureMode(params, cli_options, literature_profile);
    }
    if (cli_options.literature_mode && !params.enable_phi_correction_)
    {
        std::cout << "[ophelie] literature-mode overrides --no-phi (phi required for Jacoutot OPHELIE form)"
                  << std::endl;
        params.enable_phi_correction_ = true;
    }

    french.sigma_glass = params.sigma_glass_;
    french.frequency_hz = params.frequency_;
    french.target_joule_power = params.target_joule_power_;
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);
    logOphelieFinalParams(params, cli_options);

    if (cli_options.reload_dir.empty())
    {
        cli_options.reload_dir = "./reload";
    }
    IO::getEnvironment().resetReloadFolder(cli_options.reload_dir, true);

    StdVec<char *> filtered_argv;
    filtered_argv.reserve(filtered_arguments.size());
    for (auto &argument : filtered_arguments)
    {
        filtered_argv.push_back(const_cast<char *>(argument.c_str()));
    }
    const int filtered_ac = static_cast<int>(filtered_argv.size());
    char **filtered_av = filtered_argv.data();

    const Real dp = french.dp;
    const Real boundary_width = 3.0 * dp;
    const BoundingBoxd system_bounds = frenchReducedDomainBounds(french, boundary_width);

    const size_t relaxation_steps = cli_options.relaxation_steps > 0 ? cli_options.relaxation_steps : 400;
    const size_t relaxation_log_every =
        cli_options.relaxation_log_every > 0 ? cli_options.relaxation_log_every : 50;
    const size_t relaxation_vtp_every = cli_options.relaxation_vtp_every;

    SPHSystem sph_system(system_bounds, dp);
    sph_system.handleCommandlineOptions(filtered_ac, filtered_av);

    printFrenchReducedCaseSummary(french);
    logFrenchReducedCoilGeometry(french);

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height,
                                                                              french.glass_mesh_resolution));
    glass_body.defineAdaptation<SPHAdaptation>(1.0, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    LevelSetShape &glass_level_set = defineOphelieSolidLevelSet(glass_body);

    UniquePtr<SolidBody> coil_visual_body;
    UniquePtr<SolidBody> crucible_visual_body;
    LevelSetShape *coil_visual_level_set = nullptr;
    LevelSetShape *crucible_visual_level_set = nullptr;
    if (french.enable_coil_visual)
    {
        coil_visual_body = makeUnique<SolidBody>(
            sph_system, makeShared<OphelieFrenchCoilVisualShape>("CoilVisualBody", french));
        coil_visual_body->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        coil_visual_body->defineMatterMaterial<Solid>();
        coil_visual_level_set = &defineOphelieSolidLevelSet(*coil_visual_body);
    }
    if (french.enable_crucible_visual)
    {
        crucible_visual_body = makeUnique<SolidBody>(
            sph_system, makeShared<OphelieFrenchCrucibleWallVisualShape>("CrucibleWallVisualBody", french));
        crucible_visual_body->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        crucible_visual_body->defineMatterMaterial<Solid>();
        crucible_visual_level_set = &defineOphelieSolidLevelSet(*crucible_visual_body);
    }

    const bool use_particle_reload = !sph_system.RunParticleRelaxation() && sph_system.ReloadParticles();
    const char *particle_source = "lattice";

    if (sph_system.RunParticleRelaxation())
    {
        glass_body.generateParticles<BaseParticles, Lattice>();
        generateFrenchReducedVisualLattice(coil_visual_body.get(), crucible_visual_body.get());
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        relaxFrenchReducedBodies(sph_system, glass_body, glass_level_set, coil_visual_body.get(),
                                 coil_visual_level_set, crucible_visual_body.get(), crucible_visual_level_set,
                                 relaxation_steps, relaxation_log_every, relaxation_vtp_every);
        writeFrenchReducedReload(glass_body, coil_visual_body.get(), crucible_visual_body.get());
        std::cout << "test_3d_ophelie_french_reduced particle relaxation finished. Reload.xml (GlassBody"
                  << (coil_visual_body ? "+CoilVisualBody" : "") << (crucible_visual_body ? "+CrucibleWallVisualBody" : "")
                  << ") -> " << IO::getEnvironment().ReloadFolder() << std::endl;
        return 0;
    }

    if (use_particle_reload)
    {
        glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
        generateFrenchReducedVisualLattice(coil_visual_body.get(), crucible_visual_body.get());
        particle_source = "reload";
        std::cout << "[ophelie] loaded relaxed GlassBody from " << IO::getEnvironment().ReloadFolder()
                  << "/Reload.xml" << std::endl;
    }
    else
    {
        glass_body.generateParticles<BaseParticles, Lattice>();
        generateFrenchReducedVisualLattice(coil_visual_body.get(), crucible_visual_body.get());
        if (!cli_options.skip_relaxation)
        {
            sph_system.initializeSystemCellLinkedLists();
            sph_system.initializeSystemConfigurations();
            relaxFrenchReducedBodies(sph_system, glass_body, glass_level_set, coil_visual_body.get(),
                                     coil_visual_level_set, crucible_visual_body.get(), crucible_visual_level_set,
                                     relaxation_steps, relaxation_log_every, relaxation_vtp_every);
            particle_source = "relax";
        }
        else
        {
            std::cout << "[ophelie] --skip-relax: using lattice particles without relaxation" << std::endl;
        }
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

    OphelieFrenchEmSolveResult em_result;
    OphelieFrenchSelfInductionPicardResult picard_result;
    const bool use_picard_self_induction = params.enable_self_induction_ && params.enable_phi_correction_;
    if (use_picard_self_induction)
    {
        std::cout << "[ophelie] french_reduced: Picard self-induction (experimental; not literature_passed)"
                  << std::endl;
        picard_result = runFrenchReducedSelfInductionPicard<MainExecutionPolicy>(glass_body, *glass_inner, glass_names,
                                                                               params, french);
        em_result.phi_solver_rel_residual = picard_result.phi_solver_rel_residual;
        em_result.phi_eq_res_vol = picard_result.phi_eq_res_vol;
        em_result.joule_power_raw = picard_result.joule_power_w;
        em_result.div_j_phi =
            computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, french.glass_radius);
        em_result.div_j_level0 = em_result.div_j_phi;
    }
    else
    {
        em_result = runFrenchReducedEmPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params,
                                                                    french);
    }

    if (literature_profile.enabled && literature_profile.calibrate_coil_current)
    {
        calibrateFrenchCoilCurrentToTargetPower(french, params, em_result.joule_power_raw);
        if (use_picard_self_induction)
        {
            picard_result = runFrenchReducedSelfInductionPicard<MainExecutionPolicy>(glass_body, *glass_inner,
                                                                                   glass_names, params, french);
            em_result.phi_solver_rel_residual = picard_result.phi_solver_rel_residual;
            em_result.phi_eq_res_vol = picard_result.phi_eq_res_vol;
            em_result.joule_power_raw = picard_result.joule_power_w;
            em_result.div_j_phi = computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names,
                                                                              french.glass_radius);
            em_result.div_j_level0 = em_result.div_j_phi;
        }
        else
        {
            em_result = runFrenchReducedEmPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params,
                                                                        french);
        }
    }

    OphelieDivJMetrics div_j_level0_metrics = em_result.div_j_level0;
    OphelieDivJMetrics div_j_phi_metrics = em_result.div_j_phi;
    const Real phi_solver_rel_residual = em_result.phi_solver_rel_residual;

    const Real joule_power_raw = em_result.joule_power_raw;
    const OpheliePowerScalingFactors scaling_factors = computeOpheliePowerScalingFactors(params, joule_power_raw);
    if (params.enable_power_scaling_)
    {
        StateDynamics<MainExecutionPolicy, ScaleOphelieElectromagneticFieldsCK> scale_em_fields(
            glass_body, glass_names, scaling_factors.field_scale, scaling_factors.power_scale);
        scale_em_fields.exec();
    }

    OphelieRunMetrics metrics = collectGlassMetrics(glass_particles, glass_names);
    metrics.n_coil = 0;
    metrics.joule_power_raw = joule_power_raw;
    metrics.joule_power_scaled = hostVolWeightedSum(glass_particles, glass_names.joule_heat, n_glass);
    metrics.power_scale = scaling_factors.power_scale;
    metrics.field_scale = scaling_factors.field_scale;
    metrics.effective_current_amplitude = scaling_factors.effective_current_amplitude;
    metrics.div_j_rel_level0 = div_j_level0_metrics.div_j_rel;
    metrics.div_j_rel_phi = div_j_phi_metrics.div_j_rel > 0.0 ? div_j_phi_metrics.div_j_rel : div_j_level0_metrics.div_j_rel;
    if (metrics.div_j_rel_phi > TinyReal && metrics.div_j_rel_level0 > TinyReal)
    {
        metrics.div_j_reduction = metrics.div_j_rel_level0 / metrics.div_j_rel_phi;
    }
    const Real div_j_l2_level0 = div_j_level0_metrics.div_j_weighted_l2;
    const Real div_j_l2_phi =
        div_j_phi_metrics.div_j_weighted_l2 > 0.0 ? div_j_phi_metrics.div_j_weighted_l2 : div_j_level0_metrics.div_j_weighted_l2;
    const Real div_j_l2_reduction = div_j_l2_level0 / (div_j_l2_phi + TinyReal);
    const Real effective_ampere_turns = french.ampere_turns * scaling_factors.field_scale;
    metrics.with_phi_correction = params.enable_phi_correction_;
    metrics.phi_solver_rel_residual = phi_solver_rel_residual;
    if (use_picard_self_induction)
    {
        metrics.self_induction_j_rel_change = picard_result.final_j_rel_change;
        metrics.self_induction_iterations_used = picard_result.self_induction_iterations;
    }
    const Real phi_eq_res_vol = em_result.phi_eq_res_vol;
    if (params.enable_phi_correction_ && !use_picard_self_induction)
    {
        const std::string phi_p0_csv = cli_options.phi_p0_csv_path.empty()
                                           ? std::string("./output/ophelie_phi_p0_diagnostics.csv")
                                           : cli_options.phi_p0_csv_path;
        appendOpheliePhiP0DiagnosticCsv(phi_p0_csv, opheliePhiP0CaseLabel(params), particle_source, n_glass, dp,
                                        params, em_result.phi_p0);
    }

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(glass_body, glass_names.a_src_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.a_src_imag);
    write_states.addToWrite<Vecd>(glass_body, glass_names.b_src_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.b_src_imag);
    write_states.addToWrite<Vecd>(glass_body, glass_names.e_imag);
    write_states.addToWrite<Vecd>(glass_body, glass_names.j_imag);
    write_states.addToWrite<Real>(glass_body, glass_names.joule_heat);
    write_states.addToWrite<Real>(glass_body, glass_names.div_j_imag);
    write_states.addToWrite<Real>(glass_body, glass_names.sigma);
    if (params.enable_phi_correction_)
    {
        write_states.addToWrite<Real>(glass_body, glass_names.phi_imag);
        write_states.addToWrite<Vecd>(glass_body, glass_names.grad_phi_imag);
    }
    if (ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        write_states.addToWrite<Real>(glass_body, glass_names.joule_heat_edge);
        write_states.addToWrite<Vecd>(glass_body, glass_names.j_edge_recon_imag);
        write_states.addToWrite<Vecd>(glass_body, glass_names.e_edge_recon_imag);
        write_states.addToWrite<Vecd>(glass_body, glass_names.j_edge_recon_real);
        write_states.addToWrite<Vecd>(glass_body, glass_names.e_edge_recon_real);
        write_states.addToWrite<Real>(glass_body, glass_names.edge_flux_residual_imag);
        write_states.addToWrite<Real>(glass_body, glass_names.edge_flux_residual_real);
        write_states.addToWrite<Real>(glass_body, glass_names.joule_heat_edge_recon_complex);
    }
    glass_body.setNewlyUpdated();
    write_states.writeToFile(0);
    logOphelieOutputArtifact(IO::getEnvironment().OutputFolder() + "/GlassBody_ite_0000000000.vtp");
    if (coil_visual_body)
    {
        BodyStatesRecordingToVtp coil_vtp(*coil_visual_body);
        coil_vtp.writeToFile(0);
    }
    if (crucible_visual_body)
    {
        BodyStatesRecordingToVtp crucible_vtp(*crucible_visual_body);
        crucible_vtp.writeToFile(0);
    }

    const bool power_ok = !params.enable_power_scaling_ ||
                          (std::isfinite(metrics.joule_power_scaled) && metrics.joule_power_scaled > 0.0);
    const bool phi_residual_ok =
        !params.enable_phi_correction_ ||
        (metrics.phi_solver_rel_residual <
         10.0 * (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES
                     ? params.phi_gmres_tolerance_
                     : params.phi_solver_kind_ == OpheliePhiSolverKind::PCG ? params.phi_pcg_tolerance_
                                                                             : params.phi_jacobi_tolerance_));
    const bool demo_fields_ok = metrics.n_glass > 0 && metrics.max_a_src > 0.0 && metrics.max_b_src > 0.0 &&
                                  metrics.max_e_imag > 0.0 && metrics.max_j_imag > 0.0 && metrics.max_joule_heat > 0.0 &&
                                  metrics.min_joule_heat >= 0.0 && power_ok;
    const bool demo_passed =
        literature_profile.enabled ? demo_fields_ok : demo_fields_ok && phi_residual_ok;

    OphelieFrenchLiteratureAcceptance literature_acceptance;
    OphelieEdgeFluxLiteratureAcceptance edge_literature_acceptance;
    const bool use_edge_literature_acceptance =
        ophelieUseEdgeFluxElectromotiveRhs(params) && em_result.edge_flux_report_valid;
    const bool phi_solver_passed =
        literature_profile.enabled && params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::DivSigmaGrad &&
                params.ophelie_current_form_ != OphelieCurrentFormKind::EdgeFlux
            ? (params.enable_phi_correction_ && phi_eq_res_vol < params.phi_eq_res_vol_gate_)
            : phi_residual_ok;
    const bool divj_continuity_passed =
        use_edge_literature_acceptance ||
        !params.enable_phi_correction_ ||
        div_j_l2_reduction >= (literature_profile.enabled ? literature_profile.div_j_l2_reduction_min : 1.0);
    if (literature_profile.enabled)
    {
        if (use_edge_literature_acceptance)
        {
            edge_literature_acceptance = evaluateEdgeFluxLiteratureAcceptance(
                params, literature_profile, metrics, em_result.edge_flux_report, em_result.edge_power,
                phi_solver_rel_residual, em_result.joule_power_recon_edge, em_result.joule_power_particle,
                em_result.q_antisym, em_result.q_antisym_valid, em_result.q_spatial, em_result.q_spatial_valid);
            logEdgeFluxLiteratureAcceptance(edge_literature_acceptance);
            literature_acceptance.phi_enabled = params.enable_phi_correction_;
            literature_acceptance.power_scaling_off = !params.enable_power_scaling_;
            literature_acceptance.phi_residual_ok = edge_literature_acceptance.phi_residual_ok;
            literature_acceptance.power_raw_ok = edge_literature_acceptance.power_raw_ok;
            literature_acceptance.div_j_l2_improved = true;
            literature_acceptance.fields_ok = edge_literature_acceptance.fields_ok;
            literature_acceptance.production_operator = true;
            literature_acceptance.passed = edge_literature_acceptance.edge_stage1_residual_passed &&
                                           edge_literature_acceptance.power_raw_ok;
            literature_acceptance.production_literature_passed = edge_literature_acceptance.production_literature_passed;
        }
        else
        {
            literature_acceptance = evaluateFrenchLiteratureAcceptance(params, literature_profile, metrics,
                                                                       div_j_l2_reduction, joule_power_raw,
                                                                       phi_eq_res_vol);
            logFrenchLiteratureAcceptance(params, literature_profile, literature_acceptance, div_j_l2_reduction);
        }
    }
    const bool passed =
        literature_profile.enabled
            ? (demo_passed && phi_solver_passed &&
               (use_edge_literature_acceptance ? edge_literature_acceptance.production_literature_passed
                                               : literature_acceptance.passed && literature_acceptance.power_raw_ok))
            : demo_passed;

    std::cout << "test_3d_ophelie_french_reduced source=multiloop-filament particles=" << particle_source
              << " mode=" << (literature_profile.enabled ? "literature" : "demo")
              << " dp=" << dp << " n_glass=" << metrics.n_glass
              << " n_loops=" << french.coil.num_loops << " frequency=" << params.frequency_
              << " sigma=" << params.sigma_glass_               << " phi_correction=" << (params.enable_phi_correction_ ? 1 : 0)
              << " phi_solver_rel_res=" << metrics.phi_solver_rel_residual << " phi_eq_res_vol=" << phi_eq_res_vol;
    if (params.enable_phi_correction_ && !use_picard_self_induction)
    {
        std::cout << " rhs_integral=" << em_result.phi_p0.rhs.rhs_volume_integral
                  << " rhs_mean=" << em_result.phi_p0.rhs.rhs_volume_mean
                  << " Jn_post_phi_rel=" << em_result.phi_p0.boundary_jn_post_phi.jn_boundary_rel
                  << " phi_p0_case=" << opheliePhiP0CaseLabel(params);
    }
    std::cout << " max_ASrc=" << metrics.max_a_src
              << " max_BSrc=" << metrics.max_b_src << " max_EImag=" << metrics.max_e_imag
              << " max_JImag=" << metrics.max_j_imag << " P_raw=" << metrics.joule_power_raw
              << " P_scaled=" << metrics.joule_power_scaled << " power_scale=" << metrics.power_scale
              << " field_scale=" << metrics.field_scale << " ampere_turns_eff=" << effective_ampere_turns
              << " divJ_L0=" << metrics.div_j_rel_level0 << " divJ_phi=" << metrics.div_j_rel_phi
              << " divJ_red=" << metrics.div_j_reduction << " divJ_L2_L0=" << div_j_l2_level0
              << " divJ_L2_phi=" << div_j_l2_phi << " reconstructed_divJ_red=" << div_j_l2_reduction;
    if (use_edge_literature_acceptance)
    {
        std::cout << " edge_res_red=" << em_result.edge_flux_report.edge_res_red_l2
                  << " P_graph_edge=" << em_result.joule_power_graph_edge
                  << " P_recon_edge=" << em_result.joule_power_recon_edge
                  << " P_particle=" << em_result.joule_power_particle
                  << " P_graph_over_recon=" << em_result.edge_power.p_graph_over_recon;
    }
    if (use_picard_self_induction)
    {
        std::cout << " self_ind_iters=" << metrics.self_induction_iterations_used
                  << " self_ind_J_rel=" << metrics.self_induction_j_rel_change
                  << " A_ind_over_A_coil=" << picard_result.a_ind_over_a_coil
                  << " B_ind_over_B_coil=" << picard_result.b_ind_over_b_coil;
    }
    std::cout << " demo_passed=" << (demo_passed ? 1 : 0) << " phi_solver_passed=" << (phi_solver_passed ? 1 : 0)
              << " divJ_continuity_passed=" << (divj_continuity_passed ? 1 : 0)
              << " literature_passed=" << (literature_profile.enabled && literature_acceptance.passed ? 1 : 0)
              << " production_literature_passed="
              << (literature_profile.enabled && literature_acceptance.production_literature_passed ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
