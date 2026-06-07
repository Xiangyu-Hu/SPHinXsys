#ifndef APHI_CONTACT_GMRES_TEST_HELPERS_H
#define APHI_CONTACT_GMRES_TEST_HELPERS_H

#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.hpp"
#include "electromagnetic_dynamics/aphi_matrix_free_solve_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_bodywise_residual_helpers.h"
#include "electromagnetic_dynamics/aphi_multibody_contact_gmres_ck.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiTwoBodyInterfaceCase
{
    SPHSystem sph_system;
    SolidBody left_body;
    SolidBody right_body;
    UniquePtr<Inner<>> left_inner_ck;
    UniquePtr<Inner<>> right_inner_ck;
    UniquePtr<Contact<>> left_contact_ck;
    UniquePtr<Contact<>> right_contact_ck;

    AphiTwoBodyInterfaceCase(Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
                               int ac, char *av[])
        : sph_system(BoundingBoxd(Vecd(-boundary_width, -boundary_width, -boundary_width),
                                  Vecd(body_length + boundary_width, body_height + boundary_width,
                                       body_width + boundary_width)),
                     dp_0),
          left_body(sph_system,
                    makeShared<AphiHalfSpaceBoxShape>(
                        "LeftBody", Vecd(0.25 * body_length, 0.5 * body_height, 0.5 * body_width),
                        Vecd(0.25 * body_length, 0.5 * body_height, 0.5 * body_width))),
          right_body(sph_system,
                     makeShared<AphiHalfSpaceBoxShape>(
                         "RightBody", Vecd(0.75 * body_length, 0.5 * body_height, 0.5 * body_width),
                         Vecd(0.25 * body_length, 0.5 * body_height, 0.5 * body_width)))
    {
        if (ac > 0)
        {
            sph_system.handleCommandlineOptions(ac, av);
        }
        for (auto *body_ptr : {&left_body, &right_body})
        {
            body_ptr->defineAdaptation<SPHAdaptation>(1.15, 1.0);
            body_ptr->defineMaterial<Solid>();
            body_ptr->defineBodyLevelSetShape();
            body_ptr->generateParticles<BaseParticles, Lattice>();
        }
        left_inner_ck = makeUnique<Inner<>>(left_body);
        right_inner_ck = makeUnique<Inner<>>(right_body);
        left_contact_ck = makeUnique<Contact<>>(left_body, StdVec<RealBody *>{&right_body});
        right_contact_ck = makeUnique<Contact<>>(right_body, StdVec<RealBody *>{&left_body});
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
    }

    Inner<> &left_inner() { return *left_inner_ck; }
    Inner<> &right_inner() { return *right_inner_ck; }
    Contact<> &left_contact() { return *left_contact_ck; }
    Contact<> &right_contact() { return *right_contact_ck; }

    void updateRelations()
    {
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_left_cell_linked_list(left_body);
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_right_cell_linked_list(right_body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_left_inner(left_inner());
        UpdateRelation<MainExecutionPolicy, Inner<>> update_right_inner(right_inner());
        UpdateRelation<MainExecutionPolicy, Contact<>> update_left_contact(left_contact());
        UpdateRelation<MainExecutionPolicy, Contact<>> update_right_contact(right_contact());
        update_left_cell_linked_list.exec();
        update_right_cell_linked_list.exec();
        update_left_inner.exec();
        update_right_inner.exec();
        update_left_contact.exec();
        update_right_contact.exec();
    }
};

inline Real bodyDiscreteMmsOperatorDefect(SPHBody &body, Inner<> &inner, Contact<> &contact,
                                          const AphiVariableNames &names, const AphiLhsAssemblyOptions &options,
                                          const AphiBlockNames &approx_block, const AphiBlockNames &exact_block)
{
    StateDynamics<MainExecutionPolicy, AphiBlockLinearCombinationCK> form_difference(
        body, names.t, Real(1), Real(-1), approx_block, exact_block);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_difference(body, inner, contact, names.t, names.lhs,
                                                                         names.material, options.omega, options);
    form_difference.exec();
    apply_difference.exec();

    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Real defect_norm = hostBlockNorm(particles, names.lhs, total_real_particles);
    const Real rhs_norm_l2 = hostBlockNorm(particles, names.rhs, total_real_particles);
    return defect_norm / (rhs_norm_l2 + TinyReal);
}

inline Real twoBodyDiscreteMmsOperatorDefect(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                             const AphiLhsAssemblyOptions &options,
                                             const AphiBlockNames &approx_block, const AphiBlockNames &exact_block)
{
    const Real left_defect = bodyDiscreteMmsOperatorDefect(case_setup.left_body, case_setup.left_inner(),
                                                         case_setup.left_contact(), names, options, approx_block,
                                                         exact_block);
    const Real right_defect = bodyDiscreteMmsOperatorDefect(case_setup.right_body, case_setup.right_inner(),
                                                            case_setup.right_contact(), names, options, approx_block,
                                                            exact_block);
    return std::max(left_defect, right_defect);
}

inline Real bodyCoreRelativeBlockError(SPHBody &body, const AphiVariableNames &names, Real body_length, Real body_height,
                                     Real body_width, Real core_shell)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    return coreRelativeBlockError(particles, names.solution, names.r_hat, positions, total_real_particles, body_length,
                                body_height, body_width, core_shell);
}

inline Real twoBodyCoreRelativeBlockError(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                        Real body_length, Real body_height, Real body_width, Real core_shell)
{
    return std::max(bodyCoreRelativeBlockError(case_setup.left_body, names, body_length, body_height, body_width,
                                              core_shell),
                    bodyCoreRelativeBlockError(case_setup.right_body, names, body_length, body_height, body_width,
                                               core_shell));
}

inline Real twoBodyInterfaceBandTrueResidualNorm(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                               Real body_length, Real body_height, Real body_width, Real core_shell,
                                               Real x_interface, Real band_half_width)
{
    Real sum_squared = 0.0;
    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        const Real band_norm =
            hostInterfaceBandTrueResidualNorm(particles, names.true_residual, positions, total_real_particles,
                                              body_length, body_height, body_width, core_shell, x_interface,
                                              band_half_width);
        sum_squared += band_norm * band_norm;
    }
    return std::sqrt(sum_squared);
}

struct AphiInterfaceMmsRunMetrics
{
    AphiMatrixFreeSolverResult solver_result{};
    Real global_true_rel = 0.0;
    Real left_true_rel = 0.0;
    Real right_true_rel = 0.0;
    Real max_bodywise_true_rel = 0.0;
    Real left_rhs_norm = 0.0;
    Real right_rhs_norm = 0.0;
    Real left_solution_norm = 0.0;
    Real right_solution_norm = 0.0;
    Real left_continuous_error = 0.0;
    Real right_continuous_error = 0.0;
    Real discrete_mms_defect = 0.0;
    Real continuous_field_relative_error = 0.0;
    Real interface_band_rel = 0.0;
};

inline void setupTwoBodyInterfaceMmsFields(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                           Real sigma_conductor, Real sigma_air, Real nu, Real x_interface)
{
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(case_setup.left_body, sigma_conductor,
                                                                                nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(case_setup.right_body, sigma_air, nu,
                                                                                   names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(
        case_setup.left_body, sigma_conductor, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(case_setup.right_body, sigma_air,
                                                                                       nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_left_sigma(
        case_setup.left_body, sigma_conductor, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_right_sigma(case_setup.right_body,
                                                                                         sigma_air, names.material);
    StateDynamics<MainExecutionPolicy, benchmark::AssignInterfaceFluxMatchedPhiFieldsCK> assign_left_exact(
        case_setup.left_body, x_interface, sigma_conductor, sigma_air, names.solution);
    StateDynamics<MainExecutionPolicy, benchmark::AssignInterfaceFluxMatchedPhiFieldsCK> assign_right_exact(
        case_setup.right_body, x_interface, sigma_conductor, sigma_air, names.solution);

    initialize_left.exec();
    initialize_right.exec();
    set_left_material.exec();
    set_right_material.exec();
    assign_left_sigma.exec();
    assign_right_sigma.exec();
    assign_left_exact.exec();
    assign_right_exact.exec();
}

inline StdVec<AphiMultiBodyContactEntry> buildTwoBodyMultiBodyEntries(AphiTwoBodyInterfaceCase &case_setup)
{
    return StdVec<AphiMultiBodyContactEntry>{
        AphiMultiBodyContactEntry{case_setup.left_body, &case_setup.left_contact(), &case_setup.left_inner()},
        AphiMultiBodyContactEntry{case_setup.right_body, &case_setup.right_contact(), &case_setup.right_inner()}};
}

inline Real twoBodyGlobalTrueRelativeResidual(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names)
{
    Real sum_squared = 0.0;
    Real rhs_sum_squared = 0.0;
    for (auto *body_ptr : {&case_setup.left_body, &case_setup.right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        sum_squared += std::pow(hostBlockNorm(particles, names.true_residual, total_real_particles), 2);
        rhs_sum_squared += std::pow(hostBlockNorm(particles, names.rhs, total_real_particles), 2);
    }
    return std::sqrt(sum_squared) / (std::sqrt(rhs_sum_squared) + TinyReal);
}

inline void computeTwoBodyTrueResiduals(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                        const AphiLhsAssemblyOptions &options)
{
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_left_true_residual(
        case_setup.left_body, names.true_residual, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_right_true_residual(
        case_setup.right_body, names.true_residual, names.rhs, names.lhs);

    apply_left.exec();
    apply_right.exec();
    compute_left_true_residual.exec();
    compute_right_true_residual.exec();
}

inline void zeroGmresWorkspaceOnBody(SPHBody &body, const AphiGMRESWorkspaceNames &workspace)
{
    for (const AphiBlockNames &block : workspace.v_basis)
    {
        StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero(body, block);
        zero.exec();
    }
    for (const AphiBlockNames &block : workspace.z_basis)
    {
        StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero(body, block);
        zero.exec();
    }
}

struct AphiTwoBodyTrueRelativeBreakdown
{
    Real left_true_rel = 0.0;
    Real right_true_rel = 0.0;
    Real global_true_rel = 0.0;
    Real max_bodywise_true_rel = 0.0;
    Real left_rhs_norm = 0.0;
    Real right_rhs_norm = 0.0;
    Real left_solution_norm = 0.0;
    Real right_solution_norm = 0.0;
};

inline Real bodyTrueRelativeResidual(SPHBody &body, const AphiVariableNames &names)
{
    return bodywiseResidualEntry(body, names).true_rel;
}

inline AphiTwoBodyTrueRelativeBreakdown twoBodyTrueRelativeBreakdown(AphiTwoBodyInterfaceCase &case_setup,
                                                                     const AphiVariableNames &names)
{
    const AphiBodywiseTrueRelativeBreakdown bodywise =
        buildBodywiseTrueRelativeBreakdown({&case_setup.left_body, &case_setup.right_body}, names);

    AphiTwoBodyTrueRelativeBreakdown breakdown;
    breakdown.global_true_rel = bodywise.global_true_rel;
    breakdown.max_bodywise_true_rel = bodywise.max_bodywise_true_rel;
    breakdown.left_true_rel = bodywise.bodies[0].true_rel;
    breakdown.right_true_rel = bodywise.bodies[1].true_rel;
    breakdown.left_rhs_norm = bodywise.bodies[0].rhs_norm;
    breakdown.right_rhs_norm = bodywise.bodies[1].rhs_norm;
    breakdown.left_solution_norm = bodywise.bodies[0].solution_norm;
    breakdown.right_solution_norm = bodywise.bodies[1].solution_norm;
    return breakdown;
}

inline void polishTwoBodyContactWithBlockGsSweep(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                                 const AphiLhsAssemblyOptions &options,
                                                 const AphiMatrixFreeSolverOptions &solver_options,
                                                 const AphiGMRESWorkspaceNames &gmres_workspace,
                                                 UnsignedInt max_sweeps = 1)
{
    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> left_solver(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names, options, solver_options);
    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> right_solver(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names, options, solver_options);

    for (UnsignedInt sweep = 0; sweep < max_sweeps; ++sweep)
    {
        zeroGmresWorkspaceOnBody(case_setup.left_body, gmres_workspace);
        zeroGmresWorkspaceOnBody(case_setup.right_body, gmres_workspace);
        left_solver.solve();
        case_setup.updateRelations();

        zeroGmresWorkspaceOnBody(case_setup.left_body, gmres_workspace);
        zeroGmresWorkspaceOnBody(case_setup.right_body, gmres_workspace);
        right_solver.solve();
        case_setup.updateRelations();
    }
}

inline void fillTwoBodyInterfaceMmsMetrics(AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names,
                                           const AphiLhsAssemblyOptions &options,
                                           const AphiMatrixFreeSolverResult &solver_result, Real body_length,
                                           Real body_height, Real body_width, Real core_shell, Real x_interface,
                                           Real interface_band_half_width, AphiInterfaceMmsRunMetrics &metrics)
{
    metrics.solver_result = solver_result;
    computeTwoBodyTrueResiduals(case_setup, names, options);
    const AphiTwoBodyTrueRelativeBreakdown true_breakdown = twoBodyTrueRelativeBreakdown(case_setup, names);
    metrics.global_true_rel = true_breakdown.global_true_rel;
    metrics.left_true_rel = true_breakdown.left_true_rel;
    metrics.right_true_rel = true_breakdown.right_true_rel;
    metrics.max_bodywise_true_rel = true_breakdown.max_bodywise_true_rel;
    metrics.left_rhs_norm = true_breakdown.left_rhs_norm;
    metrics.right_rhs_norm = true_breakdown.right_rhs_norm;
    metrics.left_solution_norm = true_breakdown.left_solution_norm;
    metrics.right_solution_norm = true_breakdown.right_solution_norm;
    metrics.left_continuous_error =
        bodyCoreRelativeBlockError(case_setup.left_body, names, body_length, body_height, body_width, core_shell);
    metrics.right_continuous_error =
        bodyCoreRelativeBlockError(case_setup.right_body, names, body_length, body_height, body_width, core_shell);
    metrics.continuous_field_relative_error = std::max(metrics.left_continuous_error, metrics.right_continuous_error);
    const Real interface_band_true_residual =
        twoBodyInterfaceBandTrueResidualNorm(case_setup, names, body_length, body_height, body_width, core_shell,
                                             x_interface, interface_band_half_width);
    metrics.interface_band_rel = interface_band_true_residual / (metrics.solver_result.initial_residual_norm + TinyReal);
    metrics.discrete_mms_defect =
        twoBodyDiscreteMmsOperatorDefect(case_setup, names, options, names.solution, names.r_hat);
}

struct AphiTwoBodyPolishSweepSnapshot
{
    UnsignedInt polish_sweeps = 0;
    Real global_true_rel = 0.0;
    Real left_true_rel = 0.0;
    Real right_true_rel = 0.0;
    Real left_continuous_error = 0.0;
    Real right_continuous_error = 0.0;
};

/** Coupled multi-body Contact GMRES on flux-matched interface MMS (production path). */
inline AphiInterfaceMmsRunMetrics runTwoBodyCoupledContactInterfaceMms(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, const AphiLhsAssemblyOptions &options,
    const AphiMatrixFreeSolverOptions &solver_options, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real x_interface, Real interface_band_half_width, UnsignedInt polish_sweeps = 1)
{
    AphiInterfaceMmsRunMetrics metrics;

    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left_exact(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right_exact(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_exact_reference(case_setup.left_body, names.r_hat,
                                                                                names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_exact_reference(case_setup.right_body, names.r_hat,
                                                                                 names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                            names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_solution(case_setup.left_body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_solution(case_setup.right_body, names.solution);

    apply_left_exact.exec();
    apply_right_exact.exec();
    copy_left_exact_reference.exec();
    copy_right_exact_reference.exec();
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    zero_left_solution.exec();
    zero_right_solution.exec();
    case_setup.updateRelations();

    const UnsignedInt restart_dimension = solver_options.gmres.restart_dimension;
    const AphiGMRESWorkspaceNames gmres_workspace = buildAphiGMRESWorkspaceNames(restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_left_gmres_workspace(case_setup.left_body, restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_right_gmres_workspace(case_setup.right_body, restart_dimension);
    (void)register_left_gmres_workspace;
    (void)register_right_gmres_workspace;

    zeroGmresWorkspaceOnBody(case_setup.left_body, gmres_workspace);
    zeroGmresWorkspaceOnBody(case_setup.right_body, gmres_workspace);

    const StdVec<AphiMultiBodyContactEntry> entries = buildTwoBodyMultiBodyEntries(case_setup);
    AphiMultiBodyContactGMRESSolverCK<MainExecutionPolicy> solver(entries, names, gmres_workspace, options,
                                                                  solver_options.gmres);
    metrics.solver_result = toMatrixFreeSolverResult(solver.solve());
    case_setup.updateRelations();

    if (polish_sweeps > 0)
    {
        polishTwoBodyContactWithBlockGsSweep(case_setup, names, options, solver_options, gmres_workspace,
                                             polish_sweeps);
    }

    fillTwoBodyInterfaceMmsMetrics(case_setup, names, options, metrics.solver_result, body_length, body_height,
                                   body_width, core_shell, x_interface, interface_band_half_width, metrics);
    return metrics;
}

/** Informational: coupled GMRES then incremental block-GS polish sweep curve. */
inline StdVec<AphiTwoBodyPolishSweepSnapshot> runTwoBodyCoupledPolishSweepCurve(
    AphiTwoBodyInterfaceCase &case_setup, const AphiVariableNames &names, const AphiLhsAssemblyOptions &options,
    const AphiMatrixFreeSolverOptions &solver_options, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real x_interface, Real interface_band_half_width, UnsignedInt max_polish_sweeps,
    const AphiMatrixFreeSolverOptions *polish_solver_options = nullptr)
{
    const AphiMatrixFreeSolverOptions &polish_options =
        polish_solver_options != nullptr ? *polish_solver_options : solver_options;
    AphiInterfaceMmsRunMetrics metrics = runTwoBodyCoupledContactInterfaceMms(
        case_setup, names, options, solver_options, body_length, body_height, body_width, core_shell, x_interface,
        interface_band_half_width, 0);

    StdVec<AphiTwoBodyPolishSweepSnapshot> curve;
    curve.reserve(static_cast<size_t>(max_polish_sweeps) + 1);

    auto push_snapshot = [&](UnsignedInt sweep_count) {
        AphiTwoBodyPolishSweepSnapshot snapshot;
        snapshot.polish_sweeps = sweep_count;
        snapshot.global_true_rel = metrics.global_true_rel;
        snapshot.left_true_rel = metrics.left_true_rel;
        snapshot.right_true_rel = metrics.right_true_rel;
        snapshot.left_continuous_error = metrics.left_continuous_error;
        snapshot.right_continuous_error = metrics.right_continuous_error;
        curve.push_back(snapshot);
    };

    push_snapshot(0);

    const AphiGMRESWorkspaceNames gmres_workspace =
        buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension);
    for (UnsignedInt sweep = 1; sweep <= max_polish_sweeps; ++sweep)
    {
        polishTwoBodyContactWithBlockGsSweep(case_setup, names, options, polish_options, gmres_workspace, 1);
        fillTwoBodyInterfaceMmsMetrics(case_setup, names, options, metrics.solver_result, body_length, body_height,
                                       body_width, core_shell, x_interface, interface_band_half_width, metrics);
        push_snapshot(sweep);
    }
    return curve;
}

/** Legacy block-GS + per-body Contact GMRES (diagnostic baseline). */
inline AphiInterfaceMmsRunMetrics runTwoBodyContactInterfaceMms(AphiTwoBodyInterfaceCase &case_setup,
                                                                const AphiVariableNames &names,
                                                                const AphiLhsAssemblyOptions &options,
                                                                const AphiMatrixFreeSolverOptions &solver_options,
                                                                Real body_length, Real body_height, Real body_width,
                                                                Real core_shell, Real x_interface,
                                                                Real interface_band_half_width)
{
    AphiInterfaceMmsRunMetrics metrics;

    // Coupled RHS assembly: both bodies must keep exact neighbor fields during apply.
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left_exact(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right_exact(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_exact_reference(case_setup.left_body, names.r_hat,
                                                                                names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_exact_reference(case_setup.right_body, names.r_hat,
                                                                                 names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_left_lhs_to_rhs(case_setup.left_body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_right_lhs_to_rhs(case_setup.right_body, names.rhs,
                                                                            names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_solution(case_setup.left_body, names.solution);

    apply_left_exact.exec();
    apply_right_exact.exec();
    copy_left_exact_reference.exec();
    copy_right_exact_reference.exec();
    copy_left_lhs_to_rhs.exec();
    copy_right_lhs_to_rhs.exec();
    zero_left_solution.exec();
    syncAphiBlockToDevice(case_setup.left_body.getBaseParticles(), names.solution);
    case_setup.updateRelations();

    const UnsignedInt restart_dimension = solver_options.gmres.restart_dimension;
    const AphiGMRESWorkspaceNames gmres_workspace = buildAphiGMRESWorkspaceNames(restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_left_gmres_workspace(case_setup.left_body, restart_dimension);
    RegisterAphiGMRESWorkspaceCK register_right_gmres_workspace(case_setup.right_body, restart_dimension);
    (void)register_left_gmres_workspace;
    (void)register_right_gmres_workspace;

    zeroGmresWorkspaceOnBody(case_setup.left_body, gmres_workspace);
    syncAphiBlockToDevice(case_setup.left_body.getBaseParticles(), names.solution);

    // Block Gauss-Seidel: right body keeps exact/r_hat on device; only zero left for initial guess.
    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> left_solver(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names, options, solver_options);
    metrics.solver_result = left_solver.solve();

    // Right body still carries exact/r_hat from assembly (never zeroed). Left holds converged solution.
    case_setup.updateRelations();
    zeroGmresWorkspaceOnBody(case_setup.right_body, gmres_workspace);
    AphiMatrixFreeContactSolveCK<MainExecutionPolicy> right_solver(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names, options, solver_options);
    const AphiMatrixFreeSolverResult right_result = right_solver.solve();
    if (right_result.outer_iteration_count > metrics.solver_result.outer_iteration_count)
    {
        metrics.solver_result.outer_iteration_count = right_result.outer_iteration_count;
    }
    metrics.solver_result.arnoldi_step_count += right_result.arnoldi_step_count;

    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_left_solution(
        case_setup.left_body, case_setup.left_inner(), case_setup.left_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> apply_right_solution(
        case_setup.right_body, case_setup.right_inner(), case_setup.right_contact(), names.solution, names.lhs,
        names.material, options.omega, options);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_left_true_residual(
        case_setup.left_body, names.true_residual, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiComputeBlockResidualCK> compute_right_true_residual(
        case_setup.right_body, names.true_residual, names.rhs, names.lhs);

    apply_left_solution.exec();
    apply_right_solution.exec();
    compute_left_true_residual.exec();
    compute_right_true_residual.exec();
    case_setup.updateRelations();

    metrics.discrete_mms_defect =
        twoBodyDiscreteMmsOperatorDefect(case_setup, names, options, names.solution, names.r_hat);
    metrics.continuous_field_relative_error =
        twoBodyCoreRelativeBlockError(case_setup, names, body_length, body_height, body_width, core_shell);
    const Real interface_band_true_residual =
        twoBodyInterfaceBandTrueResidualNorm(case_setup, names, body_length, body_height, body_width, core_shell,
                                             x_interface, interface_band_half_width);
    metrics.interface_band_rel = interface_band_true_residual / (metrics.solver_result.initial_residual_norm + TinyReal);
    (void)right_result;
    return metrics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_GMRES_TEST_HELPERS_H
