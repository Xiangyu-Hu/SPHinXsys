/**
 * @file test_3d_ophelie_french_reduced.cpp
 * @brief French reduced cold-crucible case: cylindrical glass + prescribed multiloop filament coil.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_progress.h"
#include "sphinxsys.h"

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

} // namespace

int main(int ac, char *av[])
{
    logOphelieRunContext();

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);

    StdVec<std::string> french_filtered = filterFrenchReducedCommandLine(ac, av, french);
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
    french.sigma_glass = params.sigma_glass_;
    french.frequency_hz = params.frequency_;
    french.target_joule_power = params.target_joule_power_;
    syncFrenchReducedToParameters(french, params);

    StdVec<char *> filtered_argv;
    filtered_argv.reserve(filtered_arguments.size());
    for (auto &argument : filtered_arguments)
    {
        filtered_argv.push_back(const_cast<char *>(argument.c_str()));
    }
    const int filtered_ac = static_cast<int>(filtered_argv.size());
    char **filtered_av = filtered_argv.data();

    const Real dp = 0.05;
    const Real boundary_width = 3.0 * dp;
    const BoundingBoxd system_bounds = frenchReducedDomainBounds(french, boundary_width);

    SPHSystem sph_system(system_bounds, dp);
    sph_system.handleCommandlineOptions(filtered_ac, filtered_av);

    printFrenchReducedCaseSummary(french);

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();

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
    applyMultiloopFilamentBiotToGlass(glass_body.getBaseParticles(), glass_names, french.coil, params.mu0_,
                                      params.softening_length_);

    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        glass_body, glass_names);
    combine_vector_potential.exec();

    BaseParticles &glass_particles = glass_body.getBaseParticles();
    const size_t n_glass = glass_particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> compute_ejq_no_phi(glass_body, glass_names,
                                                                                              params);
    OphelieDivJMetrics div_j_level0_metrics;
    OphelieDivJMetrics div_j_phi_metrics;
    Real phi_solver_rel_residual = 0.0;

    if (params.enable_phi_correction_)
    {
        compute_ejq_no_phi.exec();
        div_j_level0_metrics = computeOphelieDivJImag<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, french.glass_radius);

        phi_solver_rel_residual = solvePhiImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
        InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
            *glass_inner, glass_names);
        StateDynamics<MainExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, glass_names,
                                                                                            params);
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        compute_grad_phi.exec();
        compute_ejq_with_phi.exec();
        div_j_phi_metrics = computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names,
                                                                        french.glass_radius);
    }
    else
    {
        compute_ejq_no_phi.exec();
        div_j_level0_metrics = computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names,
                                                                           french.glass_radius);
    }

    const Real joule_power_raw = hostVolWeightedSum(glass_particles, glass_names.joule_heat, n_glass);
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
    metrics.with_phi_correction = params.enable_phi_correction_;
    metrics.phi_solver_rel_residual = phi_solver_rel_residual;

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(glass_body, glass_names.a_src_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.b_src_real);
    write_states.addToWrite<Vecd>(glass_body, glass_names.e_imag);
    write_states.addToWrite<Vecd>(glass_body, glass_names.j_imag);
    write_states.addToWrite<Real>(glass_body, glass_names.joule_heat);
    write_states.addToWrite<Real>(glass_body, glass_names.sigma);
    if (params.enable_phi_correction_)
    {
        write_states.addToWrite<Real>(glass_body, glass_names.phi_imag);
        write_states.addToWrite<Vecd>(glass_body, glass_names.grad_phi_imag);
    }
    write_states.writeToFile(0);

    const bool power_ok = !params.enable_power_scaling_ ||
                          (std::isfinite(metrics.joule_power_scaled) && metrics.joule_power_scaled > 0.0);
    const bool phi_residual_ok =
        !params.enable_phi_correction_ ||
        (metrics.phi_solver_rel_residual <
         10.0 * (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES
                     ? params.phi_gmres_tolerance_
                     : params.phi_solver_kind_ == OpheliePhiSolverKind::PCG ? params.phi_pcg_tolerance_
                                                                             : params.phi_jacobi_tolerance_));
    const bool passed = metrics.n_glass > 0 && metrics.max_a_src > 0.0 && metrics.max_b_src > 0.0 &&
                        metrics.max_e_imag > 0.0 && metrics.max_j_imag > 0.0 && metrics.max_joule_heat > 0.0 &&
                        metrics.min_joule_heat >= 0.0 && power_ok && phi_residual_ok;

    std::cout << "test_3d_ophelie_french_reduced source=multiloop-filament dp=" << dp << " n_glass=" << metrics.n_glass
              << " n_loops=" << french.coil.num_loops << " frequency=" << params.frequency_
              << " sigma=" << params.sigma_glass_ << " phi_correction=" << (params.enable_phi_correction_ ? 1 : 0)
              << " phi_rel_res=" << metrics.phi_solver_rel_residual << " max_ASrc=" << metrics.max_a_src
              << " max_BSrc=" << metrics.max_b_src << " max_EImag=" << metrics.max_e_imag
              << " max_JImag=" << metrics.max_j_imag << " P_raw=" << metrics.joule_power_raw
              << " P_scaled=" << metrics.joule_power_scaled << " power_scale=" << metrics.power_scale
              << " field_scale=" << metrics.field_scale << " divJ_L0=" << metrics.div_j_rel_level0
              << " divJ_phi=" << metrics.div_j_rel_phi << " divJ_red=" << metrics.div_j_reduction
              << " divJ_L2_L0=" << div_j_l2_level0 << " divJ_L2_phi=" << div_j_l2_phi
              << " divJ_L2_red=" << div_j_l2_reduction << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
