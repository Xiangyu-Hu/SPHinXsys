#ifndef APHI_MULTIBODY_CONTACT_GMRES_CK_HPP
#define APHI_MULTIBODY_CONTACT_GMRES_CK_HPP

#include "electromagnetic_dynamics/aphi_multibody_contact_gmres_ck.h"
#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/aphi_krylov_diagnostics_ck.h"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.hpp"

#include <Eigen/Dense>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{

template <class ExecutionPolicy>
inline AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::AphiMultiBodyContactGMRESSolverCK(
    const StdVec<AphiMultiBodyContactEntry> &bodies, const AphiVariableNames &variable_names,
    const AphiGMRESWorkspaceNames &workspace_names, const AphiLhsAssemblyOptions &operator_options,
    const AphiGMRESSolverOptions &solver_options)
    : bodies_(bodies), names_(variable_names), workspace_(workspace_names), operator_options_(operator_options),
      solver_options_(solver_options)
{
}

template <class ExecutionPolicy>
inline Real AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::globalNormSquared(
    const AphiBlockNames &block_names) const
{
    Real sum = 0.0;
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared(entry.body, block_names);
        sum += norm_squared.exec();
    }
    if (solver_options_.multibody_inner_product == AphiMultiBodyInnerProductKind::EqualPerBodyAverage &&
        !bodies_.empty())
    {
        sum /= static_cast<Real>(bodies_.size());
    }
    return sum;
}

template <class ExecutionPolicy>
inline Real AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::globalDotProduct(const AphiBlockNames &block_x,
                                                                                 const AphiBlockNames &block_y) const
{
    Real sum = 0.0;
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_product(entry.body, block_x, block_y);
        sum += dot_product.exec();
    }
    if (solver_options_.multibody_inner_product == AphiMultiBodyInnerProductKind::EqualPerBodyAverage &&
        !bodies_.empty())
    {
        sum /= static_cast<Real>(bodies_.size());
    }
    return sum;
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::computeJacobiDiagonalAll()
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        if (entry.compute_jacobi_override)
        {
            entry.compute_jacobi_override();
        }
        else
        {
            AphiComputeBlockJacobiContactDynamicsBundle<ExecutionPolicy> compute_jacobi(
                entry.body, *entry.inner, *entry.contact, names_.material, operator_options_.omega, operator_options_);
            compute_jacobi.exec();
        }
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::applySolutionToLhsAll()
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        if (entry.apply_z_override)
        {
            entry.apply_z_override(names_.solution, names_.lhs);
        }
        else
        {
            AphiApplyContactDynamicsBundle<ExecutionPolicy> apply(
                entry.body, *entry.inner, *entry.contact, names_.solution, names_.lhs, names_.material,
                operator_options_.omega, operator_options_);
            apply.exec();
        }
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::computeRecursiveResidualAll()
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiComputeResidualCK> compute_recursive_residual(entry.body, names_);
        compute_recursive_residual.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::computeTrueResidualAll()
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiComputeBlockResidualCK> compute_true_residual(
            entry.body, names_.true_residual, names_.rhs, names_.lhs);
        compute_true_residual.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::applyFullContactZToSearchAll(
    const AphiBlockNames &z_block)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        if (entry.apply_z_override)
        {
            entry.apply_z_override(z_block, names_.search);
        }
        else
        {
            AphiApplyContactDynamicsBundle<ExecutionPolicy> apply(entry.body, *entry.inner, *entry.contact, z_block,
                                                                  names_.search, names_.material, operator_options_.omega,
                                                                  operator_options_);
            apply.exec();
        }
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::preconditionVToZAll(const AphiBlockNames &v_block,
                                                                                    const AphiBlockNames &z_block)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiApplyBlockJacobiInverseCK> precondition_v_to_z(
            entry.body, v_block, z_block, names_.material, operator_options_.omega, operator_options_,
            solver_options_.block_jacobi_kind);
        precondition_v_to_z.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::copyVToZAll(const AphiBlockNames &v_block,
                                                                             const AphiBlockNames &z_block)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_v_to_z(entry.body, z_block, v_block);
        copy_v_to_z.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::orthogonalizeSearchAll(Real scale_search, Real scale_v,
                                                                                        const AphiBlockNames &v_block)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> orthogonalize_w(
            entry.body, names_.search, scale_search, scale_v, names_.search, v_block);
        orthogonalize_w.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::scaleCopyAll(const AphiBlockNames &dst_block, Real alpha,
                                                                             const AphiBlockNames &src_block)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiBlockScaleCopyCK> scale_copy(entry.body, dst_block, alpha, src_block);
        scale_copy.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::axpySolutionAll(Real alpha,
                                                                                const AphiBlockNames &z_block)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiBlockAXPYCK> update_solution(entry.body, names_.solution, alpha, z_block);
        update_solution.exec();
    }
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::normalizeV0All(Real inv_beta)
{
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiBlockScaleCopyCK> normalize_v0(entry.body, workspace_.v_basis[0], inv_beta,
                                                                          names_.residual);
        normalize_v0.exec();
    }
}

template <class ExecutionPolicy>
inline AphiGMRESResult AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::finalizeBreakdown(AphiGMRESResult &result,
                                                                                              AphiGMRESBreakdownCode code,
                                                                                              Real residual_norm,
                                                                                              Real initial_norm)
{
    result.breakdown = true;
    result.breakdown_code = code;
    result.final_residual_norm = residual_norm;
    result.final_relative_residual = residual_norm / initial_norm;
    return result;
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::maybeRecomputeTrueResidual(
    UnsignedInt arnoldi_step, Real initial_norm, Real &true_relative_residual, Real &recursive_true_gap,
    Real &true_residual_norm)
{
    if (!solver_options_.recompute_true_residual)
    {
        return;
    }
    if (solver_options_.true_residual_check_interval == 0 ||
        arnoldi_step % solver_options_.true_residual_check_interval != 0)
    {
        return;
    }

    applySolutionToLhsAll();
    computeTrueResidualAll();
    const Real true_norm = std::sqrt(globalNormSquared(names_.true_residual));
    true_relative_residual = true_norm / initial_norm;
    true_residual_norm = true_norm;

    Real gap_sum_squared = 0.0;
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_residual_gap(
            entry.body, names_.t, Real(1), Real(-1), names_.residual, names_.true_residual);
        form_residual_gap.exec();
        ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_gap(entry.body, names_.t);
        gap_sum_squared += norm_squared_gap.exec();
    }
    const Real gap_norm = std::sqrt(gap_sum_squared);
    recursive_true_gap = gap_norm / (true_norm + TinyReal);
}

template <class ExecutionPolicy>
inline void AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::recomputeFinalTrueResidual(AphiGMRESResult &result,
                                                                                            Real initial_norm)
{
    applySolutionToLhsAll();
    computeTrueResidualAll();
    const Real true_norm = std::sqrt(globalNormSquared(names_.true_residual));
    result.final_true_residual_norm = true_norm;
    result.final_true_relative_residual = true_norm / initial_norm;

    Real gap_sum_squared = 0.0;
    for (const AphiMultiBodyContactEntry &entry : bodies_)
    {
        StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_residual_gap(
            entry.body, names_.t, Real(1), Real(-1), names_.residual, names_.true_residual);
        form_residual_gap.exec();
        ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_gap(entry.body, names_.t);
        gap_sum_squared += norm_squared_gap.exec();
    }
    const Real gap_norm = std::sqrt(gap_sum_squared);
    result.final_recursive_true_gap = gap_norm / (true_norm + TinyReal);
}

template <class ExecutionPolicy>
inline AphiGMRESResult AphiMultiBodyContactGMRESSolverCK<ExecutionPolicy>::solve()
{
    AphiGMRESResult result;
    const UnsignedInt restart_dimension = solver_options_.restart_dimension;

    if (!aphiGMRESRestartDimensionValid(restart_dimension) ||
        workspace_.v_basis.size() < restart_dimension + 1 || workspace_.z_basis.size() < restart_dimension)
    {
        result.breakdown = true;
        result.breakdown_code = AphiGMRESBreakdownCode::LeastSquaresFailed;
        return result;
    }

    if (solver_options_.use_block_jacobi_preconditioner)
    {
        computeJacobiDiagonalAll();
    }

    applySolutionToLhsAll();
    computeRecursiveResidualAll();

    Real initial_norm = std::sqrt(globalNormSquared(names_.residual));
    result.initial_residual_norm = initial_norm;
    if (!aphiIsFinite(initial_norm))
    {
        result.breakdown = true;
        result.breakdown_code = AphiGMRESBreakdownCode::InitialResidualNonFinite;
        return result;
    }
    if (initial_norm <= solver_options_.absolute_tolerance || initial_norm <= solver_options_.relative_tolerance)
    {
        result.final_residual_norm = initial_norm;
        result.final_relative_residual = Real(0);
        result.converged = true;
        recomputeFinalTrueResidual(result, initial_norm);
        return result;
    }

    UnsignedInt total_arnoldi_steps = 0;
    Real previous_outer_beta = initial_norm;

    for (UnsignedInt outer = 0; outer < solver_options_.max_outer_iterations; ++outer)
    {
        applySolutionToLhsAll();
        computeRecursiveResidualAll();

        const Real beta = std::sqrt(globalNormSquared(names_.residual));
        if (!aphiIsFinite(beta))
        {
            return finalizeBreakdown(result, AphiGMRESBreakdownCode::ResidualNonFinite, beta, initial_norm);
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
        maybeRecomputeTrueResidual(total_arnoldi_steps, initial_norm, true_relative_residual, recursive_true_gap,
                                   true_residual_norm);
        result.final_true_relative_residual = true_relative_residual;
        result.final_recursive_true_gap = recursive_true_gap;
        result.final_true_residual_norm = true_residual_norm;

        if (solver_options_.enable_debug_log)
        {
            std::cout << "[AphiMultiBodyContactGMRES] outer=" << outer << " beta=" << beta
                      << " rel_res=" << relative_residual << " true_rel_res=" << true_relative_residual << std::endl;
        }

        if (relative_residual < solver_options_.relative_tolerance || beta < solver_options_.absolute_tolerance)
        {
            result.converged = true;
            recomputeFinalTrueResidual(result, initial_norm);
            return result;
        }
        if (aphiIsResidualExplosion(relative_residual, solver_options_.residual_explosion_factor))
        {
            return finalizeBreakdown(result, AphiGMRESBreakdownCode::ResidualExplosion, beta, initial_norm);
        }

        normalizeV0All(Real(1) / (beta + TinyReal));

        Eigen::MatrixXd hessenberg =
            Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(restart_dimension + 1),
                                  static_cast<Eigen::Index>(restart_dimension));
        UnsignedInt arnoldi_steps = 0;

        for (UnsignedInt j = 0; j < restart_dimension; ++j)
        {
            const AphiBlockNames &v_j = workspace_.v_basis[j];
            const AphiBlockNames &z_j = workspace_.z_basis[j];

            if (solver_options_.use_block_jacobi_preconditioner)
            {
                preconditionVToZAll(v_j, z_j);
            }
            else
            {
                copyVToZAll(v_j, z_j);
            }

            applyFullContactZToSearchAll(z_j);

            for (UnsignedInt i = 0; i <= j; ++i)
            {
                const Real h_ij = globalDotProduct(names_.search, workspace_.v_basis[i]);
                hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = h_ij;
                orthogonalizeSearchAll(Real(1), -h_ij, workspace_.v_basis[i]);
            }

            const Real h_next_j = std::sqrt(globalNormSquared(names_.search));
            hessenberg(static_cast<Eigen::Index>(j + 1), static_cast<Eigen::Index>(j)) = h_next_j;

            if (!aphiIsFinite(h_next_j))
            {
                return finalizeBreakdown(result, AphiGMRESBreakdownCode::ArnoldiNormNonFinite, beta, initial_norm);
            }

            arnoldi_steps = j + 1;
            ++total_arnoldi_steps;
            result.arnoldi_step_count = total_arnoldi_steps;

            if (h_next_j <= solver_options_.happy_breakdown_tol)
            {
                break;
            }

            scaleCopyAll(workspace_.v_basis[j + 1], Real(1) / h_next_j, names_.search);
        }

        if (arnoldi_steps == 0)
        {
            return finalizeBreakdown(result, AphiGMRESBreakdownCode::HappyBreakdown, beta, initial_norm);
        }

        Eigen::MatrixXd hessenberg_sub =
            hessenberg.block(0, 0, static_cast<Eigen::Index>(arnoldi_steps + 1), static_cast<Eigen::Index>(arnoldi_steps));
        Eigen::VectorXd rhs_vector = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(arnoldi_steps + 1));
        rhs_vector(0) = beta;

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(hessenberg_sub);
        const Eigen::VectorXd y = qr.solve(rhs_vector);
        if (!y.allFinite())
        {
            return finalizeBreakdown(result, AphiGMRESBreakdownCode::LeastSquaresFailed, beta, initial_norm);
        }

        for (UnsignedInt k = 0; k < arnoldi_steps; ++k)
        {
            axpySolutionAll(y(static_cast<Eigen::Index>(k)), workspace_.z_basis[k]);
        }
    }

    applySolutionToLhsAll();
    computeRecursiveResidualAll();
    const Real final_norm = std::sqrt(globalNormSquared(names_.residual));
    result.final_residual_norm = final_norm;
    result.final_relative_residual = final_norm / initial_norm;
    recomputeFinalTrueResidual(result, initial_norm);

    if (result.final_relative_residual < solver_options_.relative_tolerance ||
        final_norm < solver_options_.absolute_tolerance)
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

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_MULTIBODY_CONTACT_GMRES_CK_HPP
