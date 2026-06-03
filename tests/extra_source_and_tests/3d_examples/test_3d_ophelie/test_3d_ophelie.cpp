/**
 * @file test_3d_ophelie.cpp
 * @brief OPHELIE-like induction: box glass + coil shell, Biot-Savart, optional PhiImag correction, E/J/Q, VTP.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_progress.h"
#include "simple_algorithms_ck.h"
#include "sphinxsys.h"

#include <iostream>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

class OphelieGlassBoxShape : public ComplexShape
{
  public:
    OphelieGlassBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class OphelieCoilShellBoxShape : public ComplexShape
{
  public:
    OphelieCoilShellBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &outer_halfsize,
                             const Vecd &inner_halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), outer_halfsize);
        subtract<GeometricShapeBox>(Transform(center), inner_halfsize);
    }
};

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
    OphelieTestCliOptions cli_options;
    const StdVec<std::string> filtered_arguments = filterOphelieTestCommandLine(ac, av, params, cli_options);
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

    const Vecd glass_halfsize(params.glass_halfsize_x_, params.glass_halfsize_y_, params.glass_halfsize_z_);
    const Vecd coil_outer_halfsize(params.coil_outer_halfsize_x_, params.coil_outer_halfsize_y_,
                                   params.coil_outer_halfsize_z_);
    const Vecd coil_inner_halfsize(params.coil_inner_halfsize_x_, params.coil_inner_halfsize_y_,
                                   params.coil_inner_halfsize_z_);

    const Vecd domain_min =
        params.glass_center_.cwiseMin(params.coil_center_ - coil_outer_halfsize) -
        Vecd(boundary_width, boundary_width, boundary_width);
    const Vecd domain_max =
        params.glass_center_.cwiseMax(params.coil_center_ + coil_outer_halfsize) +
        Vecd(boundary_width, boundary_width, boundary_width);
    const BoundingBoxd system_bounds(domain_min, domain_max);

    SPHSystem sph_system(system_bounds, dp);
    sph_system.handleCommandlineOptions(filtered_ac, filtered_av);

    SolidBody glass_body(sph_system, makeShared<OphelieGlassBoxShape>("GlassBody", params.glass_center_, glass_halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();

    SolidBody coil_body(sph_system,
                        makeShared<OphelieCoilShellBoxShape>("CoilSourceBody", params.coil_center_, coil_outer_halfsize,
                                                             coil_inner_halfsize));
    coil_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    coil_body.defineMatterMaterial<Solid>();
    coil_body.defineBodyLevelSetShape();
    coil_body.generateParticles<BaseParticles, Lattice>();

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieCoilFieldNames coil_names;
    OphelieGlassFieldNames glass_names;

    RegisterOphelieCoilFields register_coil_fields(coil_body, coil_names);
    RegisterOphelieGlassFields register_glass_fields(glass_body, glass_names);
    (void)register_coil_fields;
    (void)register_glass_fields;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);

    StateDynamics<MainExecutionPolicy, InitializeOphelieCoilSourceCK> initialize_coil_source(
        coil_body, coil_names, params, params.coil_center_);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_glass_sigma(glass_body, glass_names,
                                                                                       params.sigma_glass_);
    StateDynamics<MainExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot_savart(
        glass_body, coil_body, glass_names, coil_names, params);

    syncCoilSourceFieldsToDevice(coil_body.getBaseParticles(), coil_names);
    syncGlassElectromagneticFieldsToDevice(glass_body.getBaseParticles(), glass_names);

    assign_glass_sigma.exec();
    initialize_coil_source.exec();
    compute_biot_savart.exec();

    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_vector_potential(
        glass_body, glass_names);
    combine_vector_potential.exec();

    OphelieRunMetrics metrics;
    OphelieRunMetrics level0_metrics;
    metrics.with_phi_correction = params.enable_phi_correction_;
    Real phi_solver_rel_residual = 0.0;
    Real self_induction_j_rel_change = 0.0;
    size_t self_induction_iterations_used = 0;

    BaseParticles &glass_particles = glass_body.getBaseParticles();
    const size_t n_glass = glass_particles.TotalRealParticles();
    const size_t n_coil = coil_body.getBaseParticles().TotalRealParticles();

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> compute_ejq_no_phi(glass_body, glass_names, params);
    const Real div_j_characteristic_length = params.glass_halfsize_x_;
    OphelieDivJMetrics div_j_level0_metrics;
    OphelieDivJMetrics div_j_phi_metrics;

    Real level0_joule_power_raw = 0.0;
    if (cli_options.compare_level0)
    {
        compute_ejq_no_phi.exec();
        div_j_level0_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_characteristic_length);
        level0_metrics = collectGlassMetrics(glass_particles, glass_names);
        level0_joule_power_raw = hostVolWeightedSum(glass_particles, glass_names.joule_heat, n_glass);
        params.enable_phi_correction_ = true;
        metrics.with_phi_correction = true;
    }

    if (params.enable_self_induction_)
    {
        self_induction_j_rel_change = runOphelieSelfInductionWithPhiSolve<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, phi_solver_rel_residual, self_induction_iterations_used);
        div_j_phi_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_characteristic_length);
    }
    else if (params.enable_phi_correction_)
    {
        phi_solver_rel_residual = solvePhiImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

        InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
            *glass_inner, glass_names);
        StateDynamics<MainExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, glass_names, params);

        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
        compute_grad_phi.exec();
        compute_ejq_with_phi.exec();
        div_j_phi_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_characteristic_length);
    }
    else if (!cli_options.compare_level0)
    {
        compute_ejq_no_phi.exec();
        div_j_level0_metrics =
            computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_characteristic_length);
    }

    const Real joule_power_raw = hostVolWeightedSum(glass_particles, glass_names.joule_heat, n_glass);
    const OpheliePowerScalingFactors scaling_factors = computeOpheliePowerScalingFactors(params, joule_power_raw);
    StateDynamics<MainExecutionPolicy, ScaleOphelieElectromagneticFieldsCK> scale_em_fields(
        glass_body, glass_names, scaling_factors.field_scale, scaling_factors.power_scale);
    scale_em_fields.exec();

    metrics = collectGlassMetrics(glass_particles, glass_names);
    if (cli_options.compare_level0)
    {
        std::cout << "ophelie_level0_vs_phi"
                  << " d_max_EImag=" << (metrics.max_e_imag - level0_metrics.max_e_imag)
                  << " d_max_JImag=" << (metrics.max_j_imag - level0_metrics.max_j_imag)
                  << " d_P_raw=" << (joule_power_raw - level0_joule_power_raw)
                  << " L0_max_EImag=" << level0_metrics.max_e_imag << " Phi_max_EImag=" << metrics.max_e_imag
                  << " L0_P_raw=" << level0_joule_power_raw << " Phi_P_raw=" << joule_power_raw << std::endl;
    }
    metrics.n_coil = n_coil;
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
    std::cout << "[ophelie] power_scaling P_raw=" << metrics.joule_power_raw << " P_target=" << params.target_joule_power_
              << " power_scale=" << metrics.power_scale << " field_scale=" << metrics.field_scale
              << " I_eff=" << metrics.effective_current_amplitude << std::endl;
    std::cout << "[ophelie] divJ_rel_level0=" << metrics.div_j_rel_level0 << " divJ_rel_phi=" << metrics.div_j_rel_phi
              << " divJ_reduction=" << metrics.div_j_reduction << std::endl;
    metrics.with_phi_correction = params.enable_phi_correction_;
    metrics.max_a_coil = hostVecdFieldMax(glass_particles, glass_names.a_coil_real, n_glass);
    metrics.max_a_ind = hostVecdFieldMax(glass_particles, glass_names.a_ind_real, n_glass);
    if (params.enable_phi_correction_ || params.enable_self_induction_)
    {
        metrics.phi_solver_rel_residual = phi_solver_rel_residual;
        metrics.phi_solver_kind = params.phi_solver_kind_;
        metrics.self_induction_j_rel_change = self_induction_j_rel_change;
        metrics.self_induction_iterations_used = self_induction_iterations_used;
        metrics.max_phi_imag = hostScalarFieldMax(glass_particles, glass_names.phi_imag, n_glass);
        metrics.max_phi_rhs_imag = hostScalarFieldMax(glass_particles, glass_names.phi_rhs_imag, n_glass);
    }

    syncVariableToHost<Vecd>(coil_body.getBaseParticles(), coil_names.j_src_real);
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<Vecd>(coil_body, coil_names.j_src_real);
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
    if (params.enable_self_induction_)
    {
        write_states.addToWrite<Vecd>(glass_body, glass_names.a_coil_real);
        write_states.addToWrite<Vecd>(glass_body, glass_names.a_ind_real);
    }
    write_states.writeToFile(0);

    const bool passed = metrics.n_glass > 0 && metrics.n_coil > 0 && metrics.max_a_src > 0.0 && metrics.max_b_src > 0.0 &&
                        metrics.max_e_imag > 0.0 && metrics.max_j_imag > 0.0 && metrics.max_joule_heat > 0.0 &&
                        metrics.min_joule_heat >= 0.0 && std::isfinite(metrics.joule_power_scaled) &&
                        metrics.joule_power_scaled > 0.0 &&
                        (!params.enable_phi_correction_ ||
                         (metrics.max_phi_imag > 0.0 && metrics.max_phi_rhs_imag > 0.0 &&
                          metrics.phi_solver_rel_residual <
                              10.0 * (params.phi_solver_kind_ == OpheliePhiSolverKind::GMRES
                                          ? params.phi_gmres_tolerance_
                                          : params.phi_solver_kind_ == OpheliePhiSolverKind::PCG
                                                ? params.phi_pcg_tolerance_
                                                : params.phi_jacobi_tolerance_) &&
                          (!params.enable_self_induction_ ||
                           metrics.self_induction_j_rel_change < params.self_induction_j_tolerance_)));

    std::cout << "test_3d_ophelie"
              << " dp=" << dp << " n_glass=" << metrics.n_glass << " n_coil=" << metrics.n_coil
              << " frequency=" << params.frequency_ << " sigma=" << params.sigma_glass_
              << " J0=" << params.equivalentCurrentDensity() << " phi_correction=" << (params.enable_phi_correction_ ? 1 : 0)
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
              << " max_Joule=" << metrics.max_joule_heat << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
