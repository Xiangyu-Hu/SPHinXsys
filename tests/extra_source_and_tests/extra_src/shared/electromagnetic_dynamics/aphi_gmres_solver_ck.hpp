#ifndef APHI_GMRES_SOLVER_CK_HPP
#define APHI_GMRES_SOLVER_CK_HPP

#include "electromagnetic_dynamics/aphi_gmres_solver_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_gmres_residual_decomposition_helpers.h"
#include "electromagnetic_dynamics/aphi_krylov_diagnostics_ck.h"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.hpp"

#include <Eigen/Dense>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace detail
{

inline AphiGMRESResult aphiGMRESFinalizeBreakdown(AphiGMRESResult &result, AphiGMRESBreakdownCode code,
                                                  Real residual_norm, Real initial_norm)
{
    result.breakdown = true;
    result.breakdown_code = code;
    result.final_residual_norm = residual_norm;
    result.final_relative_residual = residual_norm / initial_norm;
    return result;
}

template <class ExecutionPolicy, class ApplyToLhsExec, class TrueResidualExec>
inline void aphiGMRESRecomputeFinalTrueResidual(SPHBody &body, AphiGMRESResult &result, Real initial_norm,
                                                const AphiBlockNames &true_residual_block, ApplyToLhsExec &apply_to_lhs,
                                                TrueResidualExec &compute_true_residual,
                                                StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> &form_residual_gap,
                                                ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_true_r,
                                                ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_gap)
{
    apply_to_lhs.exec();
    compute_true_residual.exec();
    const Real true_norm = std::sqrt(norm_squared_true_r.exec());
    result.final_true_residual_norm = true_norm;
    result.final_true_relative_residual = true_norm / initial_norm;

    form_residual_gap.exec();
    const Real gap_norm = std::sqrt(norm_squared_gap.exec());
    result.final_recursive_true_gap = gap_norm / (true_norm + TinyReal);

    BaseParticles &particles = body.getBaseParticles();
    const auto decomp = test::hostTrueResidualBlockDecomposition(particles, true_residual_block,
                                                                 particles.TotalRealParticles());
    result.final_res_a_total_norm = decomp.res_a_total_norm;
    result.final_res_phi_total_norm = decomp.res_phi_total_norm;
    result.final_res_a_fraction = decomp.res_a_fraction;
    result.final_res_phi_fraction = decomp.res_phi_fraction;
}

template <class ExecutionPolicy, class ApplyToLhsExec, class TrueResidualExec>
inline void aphiGMRESMaybeRecomputeTrueResidual(
    const AphiGMRESSolverOptions &solver_options, UnsignedInt arnoldi_step, Real initial_norm,
    Real &true_relative_residual, Real &recursive_true_gap, Real &true_residual_norm, ApplyToLhsExec &apply_to_lhs,
    TrueResidualExec &compute_true_residual,
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> &form_residual_gap,
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_true_r,
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_gap)
{
    if (!solver_options.recompute_true_residual)
    {
        return;
    }
    if (solver_options.true_residual_check_interval == 0 ||
        arnoldi_step % solver_options.true_residual_check_interval != 0)
    {
        return;
    }

    apply_to_lhs.exec();
    compute_true_residual.exec();
    const Real true_norm = std::sqrt(norm_squared_true_r.exec());
    true_relative_residual = true_norm / initial_norm;
    true_residual_norm = true_norm;

    form_residual_gap.exec();
    const Real gap_norm = std::sqrt(norm_squared_gap.exec());
    recursive_true_gap = gap_norm / (true_norm + TinyReal);
}

template <class ExecutionPolicy, class ComputeJacobiExec, class ApplyToLhsExec, class ApplyMatVecFn>
inline AphiGMRESResult aphiGMRESSolve(
    SPHBody &body, const AphiVariableNames &names, const AphiGMRESWorkspaceNames &workspace,
    const AphiLhsAssemblyOptions &operator_options, const AphiGMRESSolverOptions &solver_options,
    ComputeJacobiExec &compute_jacobi, ApplyToLhsExec &apply_to_lhs, ApplyMatVecFn apply_matvec,
    StateDynamics<ExecutionPolicy, AphiComputeResidualCK> &compute_recursive_residual,
    StateDynamics<ExecutionPolicy, AphiComputeBlockResidualCK> &compute_true_residual,
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> &form_residual_gap,
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_r,
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_true_r,
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> &norm_squared_gap, const char *debug_tag)
{
    AphiGMRESResult result;
    const UnsignedInt restart_dimension = solver_options.restart_dimension;

    if (!aphiGMRESRestartDimensionValid(restart_dimension) || workspace.v_basis.size() < restart_dimension + 1 ||
        workspace.z_basis.size() < restart_dimension)
    {
        result.breakdown = true;
        result.breakdown_code = AphiGMRESBreakdownCode::LeastSquaresFailed;
        return result;
    }

    if (solver_options.use_block_jacobi_preconditioner)
    {
        compute_jacobi.exec();
    }

    apply_to_lhs.exec();
    compute_recursive_residual.exec();

    Real initial_norm = std::sqrt(norm_squared_r.exec());
    result.initial_residual_norm = initial_norm;
    if (!aphiIsFinite(initial_norm))
    {
        result.breakdown = true;
        result.breakdown_code = AphiGMRESBreakdownCode::InitialResidualNonFinite;
        return result;
    }
    if (initial_norm <= solver_options.absolute_tolerance || initial_norm <= solver_options.relative_tolerance)
    {
        result.final_residual_norm = initial_norm;
        result.final_relative_residual = Real(0);
        result.converged = true;
        aphiGMRESRecomputeFinalTrueResidual<ExecutionPolicy>(body, result, initial_norm, names.true_residual,
                                                             apply_to_lhs, compute_true_residual, form_residual_gap,
                                                             norm_squared_true_r, norm_squared_gap);
        return result;
    }

    UnsignedInt total_arnoldi_steps = 0;
    Real previous_outer_beta = initial_norm;

    for (UnsignedInt outer = 0; outer < solver_options.max_outer_iterations; ++outer)
    {
        apply_to_lhs.exec();
        compute_recursive_residual.exec();

        const Real beta = std::sqrt(norm_squared_r.exec());
        if (!aphiIsFinite(beta))
        {
            return aphiGMRESFinalizeBreakdown(result, AphiGMRESBreakdownCode::ResidualNonFinite, beta, initial_norm);
        }
        const Real relative_residual = beta / initial_norm;
        if (outer > 0 && beta > previous_outer_beta + TinyReal)
        {
            result.monotonic_outer_residual = false;
        }
        previous_outer_beta = beta;
        result.outer_iteration_count = outer + 1;
        result.final_residual_norm = beta;
        result.final_relative_residual = relative_residual;

        Real true_relative_residual = relative_residual;
        Real recursive_true_gap = Real(0);
        Real true_residual_norm = result.final_true_residual_norm;
        aphiGMRESMaybeRecomputeTrueResidual<ExecutionPolicy>(
            solver_options, total_arnoldi_steps, initial_norm, true_relative_residual, recursive_true_gap,
            true_residual_norm, apply_to_lhs, compute_true_residual, form_residual_gap, norm_squared_true_r,
            norm_squared_gap);
        result.final_true_relative_residual = true_relative_residual;
        result.final_recursive_true_gap = recursive_true_gap;
        result.final_true_residual_norm = true_residual_norm;

        if (solver_options.enable_debug_log)
        {
            std::cout << debug_tag << " outer=" << outer << " beta=" << beta << " rel_res=" << relative_residual
                      << " true_rel_res=" << true_relative_residual << std::endl;
        }

        if (relative_residual < solver_options.relative_tolerance || beta < solver_options.absolute_tolerance)
        {
            result.converged = true;
            aphiGMRESRecomputeFinalTrueResidual<ExecutionPolicy>(body, result, initial_norm, names.true_residual,
                                                                   apply_to_lhs, compute_true_residual, form_residual_gap,
                                                                   norm_squared_true_r, norm_squared_gap);
            return result;
        }
        if (aphiIsResidualExplosion(relative_residual, solver_options.residual_explosion_factor))
        {
            return aphiGMRESFinalizeBreakdown(result, AphiGMRESBreakdownCode::ResidualExplosion, beta, initial_norm);
        }

        {
            StateDynamics<ExecutionPolicy, AphiBlockScaleCopyCK> normalize_v0(
                body, workspace.v_basis[0], Real(1) / (beta + TinyReal), names.residual);
            normalize_v0.exec();
        }

        Eigen::MatrixXd hessenberg =
            Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(restart_dimension + 1),
                                  static_cast<Eigen::Index>(restart_dimension));
        UnsignedInt arnoldi_steps = 0;

        for (UnsignedInt j = 0; j < restart_dimension; ++j)
        {
            const AphiBlockNames &v_j = workspace.v_basis[j];
            const AphiBlockNames &z_j = workspace.z_basis[j];

            if (solver_options.use_block_jacobi_preconditioner)
            {
                StateDynamics<ExecutionPolicy, AphiApplyBlockJacobiInverseCK> precondition_v_to_z(
                    body, v_j, z_j, names.material, operator_options.omega, operator_options,
                    solver_options.block_jacobi_kind);
                precondition_v_to_z.exec();
            }
            else
            {
                StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_v_to_z(body, z_j, v_j);
                copy_v_to_z.exec();
            }

            apply_matvec(z_j);

            for (UnsignedInt i = 0; i <= j; ++i)
            {
                ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_w_vi(body, names.search, workspace.v_basis[i]);
                const Real h_ij = dot_w_vi.exec();
                hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = h_ij;

                StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> orthogonalize_w(
                    body, names.search, Real(1), -h_ij, names.search, workspace.v_basis[i]);
                orthogonalize_w.exec();
            }

            ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_w(body, names.search);
            const Real h_next_j = std::sqrt(norm_squared_w.exec());
            hessenberg(static_cast<Eigen::Index>(j + 1), static_cast<Eigen::Index>(j)) = h_next_j;

            if (!aphiIsFinite(h_next_j))
            {
                return aphiGMRESFinalizeBreakdown(result, AphiGMRESBreakdownCode::ArnoldiNormNonFinite, beta,
                                                  initial_norm);
            }

            arnoldi_steps = j + 1;
            ++total_arnoldi_steps;
            result.arnoldi_step_count = total_arnoldi_steps;

            if (h_next_j <= solver_options.happy_breakdown_tol)
            {
                break;
            }

            StateDynamics<ExecutionPolicy, AphiBlockScaleCopyCK> normalize_next_v(
                body, workspace.v_basis[j + 1], Real(1) / h_next_j, names.search);
            normalize_next_v.exec();
        }

        if (arnoldi_steps == 0)
        {
            return aphiGMRESFinalizeBreakdown(result, AphiGMRESBreakdownCode::HappyBreakdown, beta, initial_norm);
        }

        Eigen::MatrixXd hessenberg_sub =
            hessenberg.block(0, 0, static_cast<Eigen::Index>(arnoldi_steps + 1), static_cast<Eigen::Index>(arnoldi_steps));
        Eigen::VectorXd rhs_vector = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(arnoldi_steps + 1));
        rhs_vector(0) = beta;

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(hessenberg_sub);
        const Eigen::VectorXd y = qr.solve(rhs_vector);
        if (!y.allFinite())
        {
            return aphiGMRESFinalizeBreakdown(result, AphiGMRESBreakdownCode::LeastSquaresFailed, beta, initial_norm);
        }

        for (UnsignedInt k = 0; k < arnoldi_steps; ++k)
        {
            StateDynamics<ExecutionPolicy, AphiBlockAXPYCK> update_solution(
                body, names.solution, y(static_cast<Eigen::Index>(k)), workspace.z_basis[k]);
            update_solution.exec();
        }
    }

    apply_to_lhs.exec();
    compute_recursive_residual.exec();
    const Real final_norm = std::sqrt(norm_squared_r.exec());
    result.final_residual_norm = final_norm;
    result.final_relative_residual = final_norm / initial_norm;
    aphiGMRESRecomputeFinalTrueResidual<ExecutionPolicy>(body, result, initial_norm, names.true_residual, apply_to_lhs,
                                                         compute_true_residual, form_residual_gap, norm_squared_true_r,
                                                         norm_squared_gap);

    if (result.final_relative_residual < solver_options.relative_tolerance ||
        final_norm < solver_options.absolute_tolerance)
    {
        result.converged = true;
        result.breakdown = false;
        result.breakdown_code = AphiGMRESBreakdownCode::None;
        return result;
    }

    result.converged = false;
    result.breakdown = false;
    result.breakdown_code = AphiGMRESBreakdownCode::MaxOuterIterationsReached;
    return result;
}

} // namespace detail

template <class ExecutionPolicy>
inline AphiGMRESSolverCK<ExecutionPolicy>::AphiGMRESSolverCK(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiGMRESWorkspaceNames &workspace_names, const AphiLhsAssemblyOptions &operator_options,
    const AphiGMRESSolverOptions &solver_options)
    : body_(sph_body), inner_(inner_relation), names_(variable_names), workspace_(workspace_names),
      operator_options_(operator_options), solver_options_(solver_options),
      compute_jacobi_diagonal_(DynamicsArgs(inner_relation, names_.material, operator_options_.omega, operator_options_)),
      apply_solution_to_lhs_(sph_body, inner_relation, names_.solution, names_.lhs, names_.material,
                             operator_options_.omega, operator_options_),
      compute_recursive_residual_(sph_body, names_),
      compute_true_residual_(sph_body, names_.true_residual, names_.rhs, names_.lhs),
      form_residual_gap_(sph_body, names_.t, Real(1), Real(-1), names_.residual, names_.true_residual),
      norm_squared_r_(sph_body, names_.residual),
      norm_squared_true_r_(sph_body, names_.true_residual),
      norm_squared_gap_(sph_body, names_.t)
{
}

template <class ExecutionPolicy>
inline AphiGMRESResult AphiGMRESSolverCK<ExecutionPolicy>::solve()
{
    return detail::aphiGMRESSolve<ExecutionPolicy>(
        body_, names_, workspace_, operator_options_, solver_options_, compute_jacobi_diagonal_,
        apply_solution_to_lhs_,
        [&](const AphiBlockNames &z_j) {
            AphiApplyDynamicsBundle<ExecutionPolicy> apply_z_to_w(body_, inner_, z_j, names_.search, names_.material,
                                                                  operator_options_.omega, operator_options_);
            apply_z_to_w.exec();
        },
        compute_recursive_residual_, compute_true_residual_, form_residual_gap_, norm_squared_r_, norm_squared_true_r_,
        norm_squared_gap_, "[AphiGMRES]");
}

template <class ExecutionPolicy>
inline AphiGMRESContactSolverCK<ExecutionPolicy>::AphiGMRESContactSolverCK(
    SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation, const AphiVariableNames &variable_names,
    const AphiGMRESWorkspaceNames &workspace_names, const AphiLhsAssemblyOptions &operator_options,
    const AphiGMRESSolverOptions &solver_options)
    : body_(sph_body), inner_(inner_relation), contact_(contact_relation), names_(variable_names),
      workspace_(workspace_names), operator_options_(operator_options), solver_options_(solver_options),
      compute_jacobi_diagonal_(sph_body, inner_relation, contact_relation, names_.material, operator_options_.omega,
                               operator_options_),
      apply_solution_to_lhs_(sph_body, inner_relation, contact_relation, names_.solution, names_.lhs, names_.material,
                             operator_options_.omega, operator_options_),
      compute_recursive_residual_(sph_body, names_),
      compute_true_residual_(sph_body, names_.true_residual, names_.rhs, names_.lhs),
      form_residual_gap_(sph_body, names_.t, Real(1), Real(-1), names_.residual, names_.true_residual),
      norm_squared_r_(sph_body, names_.residual), norm_squared_true_r_(sph_body, names_.true_residual),
      norm_squared_gap_(sph_body, names_.t)
{
}

template <class ExecutionPolicy>
inline AphiGMRESResult AphiGMRESContactSolverCK<ExecutionPolicy>::solve()
{
    return detail::aphiGMRESSolve<ExecutionPolicy>(
        body_, names_, workspace_, operator_options_, solver_options_, compute_jacobi_diagonal_,
        apply_solution_to_lhs_,
        [&](const AphiBlockNames &z_j) {
            AphiApplyContactBlockDiagonalDynamicsBundle<ExecutionPolicy> apply_z_to_w(
                body_, inner_, contact_, z_j, names_.search, names_.material, operator_options_.omega,
                operator_options_);
            apply_z_to_w.exec();
        },
        compute_recursive_residual_, compute_true_residual_, form_residual_gap_, norm_squared_r_, norm_squared_true_r_,
        norm_squared_gap_, "[AphiGMRESContact]");
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_SOLVER_CK_HPP
