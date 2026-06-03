#ifndef ELECTROMAGNETIC_OPHELIE_PHI_GMRES_H
#define ELECTROMAGNETIC_OPHELIE_PHI_GMRES_H

#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_progress.h"

#include <Eigen/Dense>
#include <iostream>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OpheliePhiGMRESResult
{
    UnsignedInt outer_iteration_count = 0;
    UnsignedInt arnoldi_step_count = 0;
    Real initial_relative_residual = 0.0;
    Real final_relative_residual = 0.0;
    bool converged = false;
};

template <class ExecutionPolicy>
inline Real solvePhiImagGMRES(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                              const OphelieParameters &params)
{
    setupOpheliePhiImagRhsProblem<ExecutionPolicy>(glass_body, inner, names, params);

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_diag(glass_body, names.phi_laplace_diag);
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceDiagonalCK<Inner<>>> compute_diag(
        inner, names.sigma, names.phi_laplace_diag, params.pair_weight_regularization_);
    zero_diag.exec();
    compute_diag.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    const UnsignedInt restart_dimension = std::max(params.phi_gmres_restart_dimension_, UnsignedInt(1));
    const UnsignedInt max_outer_iterations = std::max(params.phi_gmres_max_outer_iterations_, UnsignedInt(1));

    StdVec<Real> solution(n, Real(0));
    StdVec<Real> rhs(n, Real(0));
    StdVec<Real> diagonal(n, Real(0));
    StdVec<Real> operator_output(n, Real(0));
    StdVec<Real> residual(n, Real(0));
    StdVec<Real> workspace(n, Real(0));
    StdVec<Real> preconditioned(n, Real(0));
    StdVec<StdVec<Real>> krylov_basis(restart_dimension + 1, StdVec<Real>(n, Real(0)));

    hostReadScalarField(particles, names.phi_rhs_imag, rhs.data(), n);
    hostReadScalarField(particles, names.phi_laplace_diag, diagonal.data(), n);

    for (size_t i = 0; i < n; ++i)
    {
        const Real preconditioner = diagonal[i] + params.phi_gauge_penalty_ + TinyReal;
        solution[i] = rhs[i] / preconditioner;
    }
    hostAssignScalarField(particles, names.phi_imag, solution.data(), n);

    auto apply_operator = [&](const Real *input, Real *output)
    {
        hostAssignScalarField(particles, names.phi_imag, input, n);
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        hostReadScalarField(particles, names.phi_lhs_imag, output, n);
    };

    auto apply_preconditioner = [&](const Real *input, Real *output)
    {
        for (size_t i = 0; i < n; ++i)
        {
            output[i] = input[i] / (diagonal[i] + params.phi_gauge_penalty_ + TinyReal);
        }
    };

    apply_operator(solution.data(), operator_output.data());
    for (size_t i = 0; i < n; ++i)
    {
        residual[i] = rhs[i] - operator_output[i];
    }

    const Real rhs_norm = hostVolWeightedNorm(particles, rhs.data(), n);
    const Real rhs_max = hostScalarFieldMax(particles, names.phi_rhs_imag, n);
    Real beta = hostVolWeightedNorm(particles, residual.data(), n);
    Real relative_residual = beta / (rhs_norm + TinyReal);

    if (!std::isfinite(relative_residual) ||
        (rhs_max > TinyReal && relative_residual < params.phi_gmres_tolerance_))
    {
        hostAssignScalarField(particles, names.phi_imag, solution.data(), n);
        return relative_residual;
    }

    StdVec<Real> hessenberg_cos(restart_dimension, Real(0));
    StdVec<Real> hessenberg_sin(restart_dimension, Real(0));

    for (UnsignedInt outer = 0; outer < max_outer_iterations; ++outer)
    {
        const Real inv_beta = Real(1.0) / (beta + TinyReal);
        for (size_t i = 0; i < n; ++i)
        {
            krylov_basis[0][i] = residual[i] * inv_beta;
        }

        Eigen::VectorXd givens_rhs = Eigen::VectorXd::Zero(restart_dimension + 1);
        givens_rhs(0) = beta;
        Eigen::MatrixXd hessenberg =
            Eigen::MatrixXd::Zero(restart_dimension + 1, restart_dimension);

        UnsignedInt krylov_dimension = 0;
        for (UnsignedInt inner_step = 0; inner_step < restart_dimension; ++inner_step)
        {
            apply_preconditioner(krylov_basis[inner_step].data(), preconditioned.data());
            apply_operator(preconditioned.data(), workspace.data());

            for (UnsignedInt i = 0; i <= inner_step; ++i)
            {
                const Real projection =
                    hostVolWeightedDot(particles, krylov_basis[i].data(), workspace.data(), n);
                hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(inner_step)) = projection;
                hostSubtractScaledVector(workspace.data(), krylov_basis[i].data(), projection, n);
            }

            const Real subdiagonal_norm = hostVolWeightedNorm(particles, workspace.data(), n);
            hessenberg(static_cast<Eigen::Index>(inner_step + 1), static_cast<Eigen::Index>(inner_step)) =
                subdiagonal_norm;

            if (!std::isfinite(subdiagonal_norm))
            {
                krylov_dimension = inner_step;
                break;
            }

            if (subdiagonal_norm < TinyReal)
            {
                krylov_dimension = inner_step + 1;
                break;
            }

            const Real inv_subdiagonal_norm = Real(1.0) / (subdiagonal_norm + TinyReal);
            for (size_t k = 0; k < n; ++k)
            {
                krylov_basis[inner_step + 1][k] = workspace[k] * inv_subdiagonal_norm;
            }

            for (UnsignedInt i = 0; i < inner_step; ++i)
            {
                const Real h0 = hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(inner_step));
                const Real h1 = hessenberg(static_cast<Eigen::Index>(i + 1), static_cast<Eigen::Index>(inner_step));
                const Real temp = hessenberg_cos[i] * h0 + hessenberg_sin[i] * h1;
                hessenberg(static_cast<Eigen::Index>(i + 1), static_cast<Eigen::Index>(inner_step)) =
                    -hessenberg_sin[i] * h0 + hessenberg_cos[i] * h1;
                hessenberg(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(inner_step)) = temp;
            }

            const Real h0 = hessenberg(static_cast<Eigen::Index>(inner_step), static_cast<Eigen::Index>(inner_step));
            const Real h1 = hessenberg(static_cast<Eigen::Index>(inner_step + 1), static_cast<Eigen::Index>(inner_step));
            const Real denom = std::sqrt(h0 * h0 + h1 * h1 + TinyReal);
            hessenberg_cos[inner_step] = h0 / denom;
            hessenberg_sin[inner_step] = h1 / denom;
            hessenberg(static_cast<Eigen::Index>(inner_step), static_cast<Eigen::Index>(inner_step)) = denom;
            hessenberg(static_cast<Eigen::Index>(inner_step + 1), static_cast<Eigen::Index>(inner_step)) = Real(0);

            const Real temp_rhs = hessenberg_cos[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step)) +
                                  hessenberg_sin[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step + 1));
            givens_rhs(static_cast<Eigen::Index>(inner_step + 1)) =
                -hessenberg_sin[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step)) +
                hessenberg_cos[inner_step] * givens_rhs(static_cast<Eigen::Index>(inner_step + 1));
            givens_rhs(static_cast<Eigen::Index>(inner_step)) = temp_rhs;

            krylov_dimension = inner_step + 1;
            relative_residual = std::abs(givens_rhs(static_cast<Eigen::Index>(inner_step + 1))) / (rhs_norm + TinyReal);
            if (relative_residual < params.phi_gmres_tolerance_)
            {
                break;
            }
        }

        if (krylov_dimension == 0)
        {
            break;
        }

        const Eigen::Index system_size = static_cast<Eigen::Index>(krylov_dimension);
        Eigen::MatrixXd hessenberg_system = hessenberg.topLeftCorner(system_size, system_size);
        Eigen::VectorXd rhs_vector = givens_rhs.head(system_size);
        Eigen::VectorXd coefficients = hessenberg_system.colPivHouseholderQr().solve(rhs_vector);

        for (Eigen::Index column = 0; column < system_size; ++column)
        {
            apply_preconditioner(krylov_basis[static_cast<size_t>(column)].data(), preconditioned.data());
            hostScaledAdd(solution.data(), coefficients(column), preconditioned.data(), n);
        }

        apply_operator(solution.data(), operator_output.data());
        for (size_t i = 0; i < n; ++i)
        {
            residual[i] = rhs[i] - operator_output[i];
        }
        beta = hostVolWeightedNorm(particles, residual.data(), n);
        relative_residual = beta / (rhs_norm + TinyReal);
        if (relative_residual < params.phi_gmres_tolerance_)
        {
            break;
        }
    }

    hostAssignScalarField(particles, names.phi_imag, solution.data(), n);
    return evaluateOpheliePhiImagRelativeResidual<ExecutionPolicy>(glass_body, inner, names, params);
}

template <class ExecutionPolicy>
inline Real solvePhiImagPCG(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                            const OphelieParameters &params)
{
    setupOpheliePhiImagRhsProblem<ExecutionPolicy>(glass_body, inner, names, params);

    StateDynamics<ExecutionPolicy, ZeroOphelieScalarFieldCK> zero_diag(glass_body, names.phi_laplace_diag);
    InteractionDynamicsCK<ExecutionPolicy, OpheliePairwiseLaplaceDiagonalCK<Inner<>>> compute_diag(
        inner, names.sigma, names.phi_laplace_diag, params.pair_weight_regularization_);
    zero_diag.exec();
    compute_diag.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StdVec<Real> solution(n, Real(0));
    StdVec<Real> rhs(n, Real(0));
    StdVec<Real> diagonal(n, Real(0));
    StdVec<Real> residual(n, Real(0));
    StdVec<Real> preconditioned_residual(n, Real(0));
    StdVec<Real> search_direction(n, Real(0));
    StdVec<Real> operator_on_search(n, Real(0));

    hostReadScalarField(particles, names.phi_rhs_imag, rhs.data(), n);
    hostReadScalarField(particles, names.phi_laplace_diag, diagonal.data(), n);

    auto apply_operator = [&](const Real *input, Real *output)
    {
        hostAssignScalarField(particles, names.phi_imag, input, n);
        applyOpheliePhiImagLhsOperator<ExecutionPolicy>(glass_body, inner, names, params);
        hostReadScalarField(particles, names.phi_lhs_imag, output, n);
    };

    auto apply_preconditioner = [&](const Real *input, Real *output)
    {
        for (size_t i = 0; i < n; ++i)
        {
            output[i] = input[i] / (diagonal[i] + params.phi_gauge_penalty_ + TinyReal);
        }
    };

    for (size_t i = 0; i < n; ++i)
    {
        solution[i] = rhs[i] / (diagonal[i] + params.phi_gauge_penalty_ + TinyReal);
    }

    apply_operator(solution.data(), operator_on_search.data());
    for (size_t i = 0; i < n; ++i)
    {
        residual[i] = rhs[i] - operator_on_search[i];
    }

    const Real rhs_norm = hostVolWeightedNorm(particles, rhs.data(), n);
    const Real rhs_max = hostScalarFieldMax(particles, names.phi_rhs_imag, n);
    Real relative_residual_l2_vol = hostVolWeightedNorm(particles, residual.data(), n) / (rhs_norm + TinyReal);
    hostAssignScalarField(particles, names.phi_imag, solution.data(), n);
    Real relative_residual_linf =
        evaluateOpheliePhiImagRelativeResidual<ExecutionPolicy>(glass_body, inner, names, params);
    if (rhs_max > TinyReal && relative_residual_linf < params.phi_pcg_tolerance_)
    {
        return relative_residual_linf;
    }

    apply_preconditioner(residual.data(), preconditioned_residual.data());
    search_direction = preconditioned_residual;
    Real rz_old = hostVolWeightedDot(particles, residual.data(), preconditioned_residual.data(), n);

    OphelieProgressLogger pcg_progress("phi_pcg");
    pcg_progress.log("n=" + std::to_string(n) + " max_iter=" + std::to_string(params.phi_pcg_max_iterations_));

    for (size_t iter = 0; iter < params.phi_pcg_max_iterations_; ++iter)
    {
        apply_operator(search_direction.data(), operator_on_search.data());
        const Real search_energy =
            hostVolWeightedDot(particles, search_direction.data(), operator_on_search.data(), n);
        if (std::abs(search_energy) < TinyReal)
        {
            break;
        }

        const Real step_length = rz_old / (search_energy + TinyReal);
        hostScaledAdd(solution.data(), step_length, search_direction.data(), n);
        for (size_t i = 0; i < n; ++i)
        {
            residual[i] -= step_length * operator_on_search[i];
        }

        relative_residual_l2_vol = hostVolWeightedNorm(particles, residual.data(), n) / (rhs_norm + TinyReal);
        if ((iter + 1) % 5 == 0 || iter + 1 == params.phi_pcg_max_iterations_ ||
            (iter + 1) % 100 == 0)
        {
            hostAssignScalarField(particles, names.phi_imag, solution.data(), n);
            relative_residual_linf =
                evaluateOpheliePhiImagRelativeResidual<ExecutionPolicy>(glass_body, inner, names, params);
            if ((iter + 1) % 100 == 0 || iter + 1 == params.phi_pcg_max_iterations_)
            {
                pcg_progress.log("iter " + std::to_string(iter + 1) + " rel_res_l2_vol=" +
                                 std::to_string(relative_residual_l2_vol) + " rel_res_linf=" +
                                 std::to_string(relative_residual_linf));
            }
            if (relative_residual_linf < params.phi_pcg_tolerance_)
            {
                break;
            }
        }

        apply_preconditioner(residual.data(), preconditioned_residual.data());
        const Real rz_new = hostVolWeightedDot(particles, residual.data(), preconditioned_residual.data(), n);
        const Real beta = rz_new / (rz_old + TinyReal);
        for (size_t i = 0; i < n; ++i)
        {
            search_direction[i] = preconditioned_residual[i] + beta * search_direction[i];
        }
        rz_old = rz_new;
    }

    hostAssignScalarField(particles, names.phi_imag, solution.data(), n);
    const Real final_residual_linf =
        evaluateOpheliePhiImagRelativeResidual<ExecutionPolicy>(glass_body, inner, names, params);
    const Real final_residual_l2_vol = hostVolWeightedNorm(particles, residual.data(), n) / (rhs_norm + TinyReal);
    pcg_progress.finish("rel_res_linf=" + std::to_string(final_residual_linf) + " rel_res_l2_vol=" +
                        std::to_string(final_residual_l2_vol));
    return final_residual_linf;
}

template <class ExecutionPolicy>
inline Real solvePhiImag(SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
                         const OphelieParameters &params)
{
    switch (params.phi_solver_kind_)
    {
    case OpheliePhiSolverKind::GMRES:
        return solvePhiImagGMRES<ExecutionPolicy>(glass_body, inner, names, params);
    case OpheliePhiSolverKind::PCG:
        return solvePhiImagPCG<ExecutionPolicy>(glass_body, inner, names, params);
    default:
        return solvePhiImagJacobi<ExecutionPolicy>(glass_body, inner, names, params);
    }
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_GMRES_H
