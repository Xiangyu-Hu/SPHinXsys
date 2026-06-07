#ifndef APHI_BICGSTAB_SOLVER_CK_HPP
#define APHI_BICGSTAB_SOLVER_CK_HPP

#include "electromagnetic_dynamics/alternate_krylov/aphi_bicgstab_solver_ck.h"
#include "electromagnetic_dynamics/aphi_krylov_diagnostics_ck.h"

#include <iostream>

namespace SPH
{
namespace electromagnetics
{

template <class ExecutionPolicy>
inline AphiBiCGStabSolverCK<ExecutionPolicy>::AphiBiCGStabSolverCK(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &operator_options, const AphiBiCGStabSolverOptions &solver_options)
    : body_(sph_body), names_(variable_names), operator_options_(operator_options), solver_options_(solver_options),
      compute_jacobi_diagonal_(DynamicsArgs(inner_relation, names_.material, operator_options_.omega, operator_options_)),
      apply_solution_to_lhs_(sph_body, inner_relation, names_.solution, names_.lhs, names_.material,
                             operator_options_.omega, operator_options_),
      apply_search_to_v_(sph_body, inner_relation, names_.search, names_.v, names_.material, operator_options_.omega,
                         operator_options_),
      apply_z_to_v_(sph_body, inner_relation, names_.z, names_.v, names_.material, operator_options_.omega,
                    operator_options_),
      apply_y_to_v_(sph_body, inner_relation, names_.y, names_.v, names_.material, operator_options_.omega,
                    operator_options_),
      apply_s_to_t_(sph_body, inner_relation, names_.s, names_.t, names_.material, operator_options_.omega,
                    operator_options_),
      compute_recursive_residual_(sph_body, names_),
      compute_true_residual_(sph_body, names_.true_residual, names_.rhs, names_.lhs),
      copy_residual_to_r_hat_(sph_body, names_.r_hat, names_.residual),
      copy_residual_to_search_(sph_body, names_.search, names_.residual),
      copy_v_to_v_old_(sph_body, names_.v_old, names_.v),
      precondition_search_to_z_(sph_body, names_.search, names_.z, names_.material, operator_options_.omega,
                                operator_options_),
      precondition_s_to_y_(sph_body, names_.s, names_.y, names_.material, operator_options_.omega, operator_options_),
      form_s_(sph_body, names_.s, Real(1), Real(0), names_.residual, names_.v),
      form_residual_gap_(sph_body, names_.t, Real(1), Real(-1), names_.residual, names_.true_residual),
      update_solution_(sph_body, names_.solution, Real(0), Real(0), names_.search, names_.s),
      update_solution_preconditioned_(sph_body, names_.solution, Real(0), Real(0), names_.z, names_.y),
      update_residual_(sph_body, names_.residual, Real(0), names_.s, names_.t),
      update_residual_preconditioned_(sph_body, names_.residual, Real(0), names_.s, names_.v),
      update_search_(sph_body, names_.search, Real(0), Real(0), names_.residual, names_.v),
      update_search_preconditioned_(sph_body, names_.search, Real(0), Real(0), names_.residual, names_.v_old),
      dot_r_hat_r_(sph_body, names_.r_hat, names_.residual),
      dot_r_hat_v_(sph_body, names_.r_hat, names_.v),
      dot_v_s_(sph_body, names_.v, names_.s),
      dot_t_s_(sph_body, names_.t, names_.s),
      dot_v_v_(sph_body, names_.v, names_.v),
      dot_t_t_(sph_body, names_.t, names_.t),
      norm_squared_r_(sph_body, names_.residual),
      norm_squared_true_r_(sph_body, names_.true_residual),
      norm_squared_gap_(sph_body, names_.t)
{
}

template <class ExecutionPolicy>
inline AphiBiCGStabSolverCK<ExecutionPolicy>::AphiBiCGStabSolverCK(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &operator_options, Real tolerance, UnsignedInt max_iterations,
    bool use_block_jacobi_preconditioner, bool enable_debug_log, Real residual_explosion_factor)
    : AphiBiCGStabSolverCK(sph_body, inner_relation, variable_names, operator_options,
                           AphiBiCGStabSolverOptions{.max_iterations = max_iterations,
                                                     .relative_tolerance = tolerance,
                                                     .use_block_jacobi_preconditioner = use_block_jacobi_preconditioner,
                                                     .enable_debug_log = enable_debug_log,
                                                     .residual_explosion_factor = residual_explosion_factor})
{
}

template <class ExecutionPolicy>
inline AphiBiCGStabResult AphiBiCGStabSolverCK<ExecutionPolicy>::finalizeBreakdown(
    AphiBiCGStabResult &result, AphiBiCGStabBreakdownCode code, Real residual_norm, Real initial_norm)
{
    result.breakdown = true;
    result.breakdown_code = code;
    result.final_residual_norm = residual_norm;
    result.final_relative_residual = residual_norm / initial_norm;
    return result;
}

template <class ExecutionPolicy>
inline void AphiBiCGStabSolverCK<ExecutionPolicy>::maybeRecomputeTrueResidual(
    UnsignedInt iteration, Real initial_norm, Real recursive_norm, Real &true_relative_residual,
    Real &recursive_true_gap, Real &true_residual_norm)
{
    if (!solver_options_.recompute_true_residual)
    {
        return;
    }
    if (solver_options_.true_residual_check_interval == 0 ||
        (iteration + 1) % solver_options_.true_residual_check_interval != 0)
    {
        return;
    }

    apply_solution_to_lhs_.exec();
    compute_true_residual_.exec();
    const Real true_norm = std::sqrt(norm_squared_true_r_.exec());
    true_relative_residual = true_norm / initial_norm;
    true_residual_norm = true_norm;

    form_residual_gap_.exec();
    const Real gap_norm = std::sqrt(norm_squared_gap_.exec());
    recursive_true_gap = gap_norm / (true_norm + TinyReal);
    (void)recursive_norm;
}

template <class ExecutionPolicy>
inline AphiBiCGStabResult AphiBiCGStabSolverCK<ExecutionPolicy>::solve()
{
    AphiBiCGStabResult result;

    if (solver_options_.use_block_jacobi_preconditioner)
    {
        compute_jacobi_diagonal_.exec();
    }

    apply_solution_to_lhs_.exec();
    compute_recursive_residual_.exec();
    copy_residual_to_r_hat_.exec();
    copy_residual_to_search_.exec();

    const Real initial_norm = std::sqrt(norm_squared_r_.exec());
    result.initial_residual_norm = initial_norm;
    if (!aphiIsFinite(initial_norm))
    {
        result.breakdown = true;
        result.breakdown_code = AphiBiCGStabBreakdownCode::InitialResidualNonFinite;
        return result;
    }
    if (initial_norm <= solver_options_.absolute_tolerance ||
        initial_norm <= solver_options_.relative_tolerance)
    {
        result.final_residual_norm = initial_norm;
        result.final_relative_residual = Real(0);
        result.converged = true;
        return result;
    }

    Real rho_old = dot_r_hat_r_.exec();
    if (!aphiIsFinite(rho_old) || aphiIsNearZeroAbs(rho_old, solver_options_.coefficient_breakdown_tol))
    {
        return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::RhoOldNearZero, initial_norm, initial_norm);
    }

    for (UnsignedInt iteration = 0; iteration < solver_options_.max_iterations; ++iteration)
    {
        if (solver_options_.use_block_jacobi_preconditioner)
        {
            precondition_search_to_z_.exec();
            apply_z_to_v_.exec();
            copy_v_to_v_old_.exec();
        }
        else
        {
            apply_search_to_v_.exec();
        }

        const Real r_hat_dot_v = dot_r_hat_v_.exec();
        if (!aphiIsFinite(r_hat_dot_v) ||
            aphiIsNearZeroAbs(r_hat_dot_v, solver_options_.coefficient_breakdown_tol))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::RHatDotVNearZero, residual_norm, initial_norm);
        }

        const Real alpha = rho_old / r_hat_dot_v;
        if (!aphiIsFinite(alpha))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::AlphaNonFinite, residual_norm, initial_norm);
        }

        form_s_.setCoefficients(Real(1), -alpha);
        form_s_.exec();

        if (solver_options_.use_block_jacobi_preconditioner)
        {
            precondition_s_to_y_.exec();
            apply_y_to_v_.exec();
        }
        else
        {
            apply_s_to_t_.exec();
        }

        const Real t_dot_s =
            solver_options_.use_block_jacobi_preconditioner ? dot_v_s_.exec() : dot_t_s_.exec();
        const Real t_dot_t =
            solver_options_.use_block_jacobi_preconditioner ? dot_v_v_.exec() : dot_t_t_.exec();
        if (!aphiIsFinite(t_dot_t) || aphiIsNearZeroAbs(t_dot_t, solver_options_.coefficient_breakdown_tol))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::TDotTNearZero, residual_norm, initial_norm);
        }

        const Real omega = t_dot_s / t_dot_t;
        if (!aphiIsFinite(omega))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::OmegaNonFinite, residual_norm, initial_norm);
        }
        if (aphiIsNearZeroAbs(omega, solver_options_.coefficient_breakdown_tol))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::OmegaNearZero, residual_norm, initial_norm);
        }

        if (solver_options_.use_block_jacobi_preconditioner)
        {
            update_solution_preconditioned_.setScalars(alpha, omega);
            update_solution_preconditioned_.exec();
            update_residual_preconditioned_.setOmega(omega);
            update_residual_preconditioned_.exec();
        }
        else
        {
            update_solution_.setScalars(alpha, omega);
            update_solution_.exec();
            update_residual_.setOmega(omega);
            update_residual_.exec();
        }

        const Real residual_norm = std::sqrt(norm_squared_r_.exec());
        const Real relative_residual = residual_norm / initial_norm;
        result.iteration_count = iteration + 1;
        result.final_residual_norm = residual_norm;
        result.final_relative_residual = relative_residual;

        Real true_relative_residual = relative_residual;
        Real recursive_true_gap = Real(0);
        Real true_residual_norm = result.final_true_residual_norm;
        maybeRecomputeTrueResidual(iteration, initial_norm, residual_norm, true_relative_residual, recursive_true_gap,
                                   true_residual_norm);
        result.final_true_relative_residual = true_relative_residual;
        result.final_recursive_true_gap = recursive_true_gap;
        result.final_true_residual_norm = true_residual_norm;

        Real beta = Real(0);
        const Real rho_new = dot_r_hat_r_.exec();
        if (iteration + 1 < solver_options_.max_iterations)
        {
            if (!aphiIsFinite(rho_new) ||
                aphiIsNearZeroAbs(rho_new, solver_options_.coefficient_breakdown_tol))
            {
                return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::RhoNewNearZero, residual_norm, initial_norm);
            }

            beta = (rho_new / rho_old) * (alpha / omega);
            if (!aphiIsFinite(beta))
            {
                return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::BetaNonFinite, residual_norm, initial_norm);
            }

            if (solver_options_.use_block_jacobi_preconditioner)
            {
                update_search_preconditioned_.setScalars(beta, omega);
                update_search_preconditioned_.exec();
            }
            else
            {
                update_search_.setScalars(beta, omega);
                update_search_.exec();
            }
            rho_old = rho_new;
        }

        if (solver_options_.enable_debug_log)
        {
            std::cout << "[AphiBiCGStab] iter=" << iteration << " rho=" << rho_new << " r_hat_dot_v=" << r_hat_dot_v
                      << " alpha=" << alpha << " t_dot_s=" << t_dot_s << " t_dot_t=" << t_dot_t << " omega=" << omega
                      << " beta=" << beta << " recursive_rel_res=" << relative_residual
                      << " true_rel_res=" << true_relative_residual << " residual_gap=" << recursive_true_gap
                      << " breakdown=" << static_cast<int>(AphiBiCGStabBreakdownCode::None) << std::endl;
        }

        if (!aphiIsFinite(residual_norm))
        {
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::ResidualNonFinite, residual_norm, initial_norm);
        }
        if (aphiIsResidualExplosion(relative_residual, solver_options_.residual_explosion_factor))
        {
            return finalizeBreakdown(result, AphiBiCGStabBreakdownCode::ResidualExplosion, residual_norm, initial_norm);
        }

        if (relative_residual < solver_options_.relative_tolerance ||
            residual_norm < solver_options_.absolute_tolerance)
        {
            result.converged = true;
            return result;
        }
    }

    result.breakdown = true;
    result.breakdown_code = AphiBiCGStabBreakdownCode::MaxIterationsReached;
    return result;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BICGSTAB_SOLVER_CK_HPP
