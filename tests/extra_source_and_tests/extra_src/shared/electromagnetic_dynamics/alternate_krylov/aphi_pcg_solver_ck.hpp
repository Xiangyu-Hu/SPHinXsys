#ifndef APHI_PCG_SOLVER_CK_HPP
#define APHI_PCG_SOLVER_CK_HPP

#include "electromagnetic_dynamics/alternate_krylov/aphi_pcg_solver_ck.h"
#include "electromagnetic_dynamics/aphi_krylov_diagnostics_ck.h"

#include <iostream>

namespace SPH
{
namespace electromagnetics
{

template <class ExecutionPolicy>
inline AphiPCGSolverCK<ExecutionPolicy>::AphiPCGSolverCK(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &operator_options, const AphiPCGSolverOptions &solver_options)
    : body_(sph_body), names_(variable_names), operator_options_(operator_options), solver_options_(solver_options),
      compute_jacobi_diagonal_(DynamicsArgs(inner_relation, names_.material, operator_options_.omega, operator_options_)),
      apply_solution_to_lhs_(sph_body, inner_relation, names_.solution, names_.lhs, names_.material,
                             operator_options_.omega, operator_options_),
      apply_search_to_v_(sph_body, inner_relation, names_.search, names_.v, names_.material, operator_options_.omega,
                         operator_options_),
      compute_recursive_residual_(sph_body, names_),
      compute_true_residual_(sph_body, names_.true_residual, names_.rhs, names_.lhs),
      copy_z_to_search_(sph_body, names_.search, names_.z),
      copy_r_to_z_(sph_body, names_.z, names_.residual),
      precondition_residual_to_z_(sph_body, names_.residual, names_.z, names_.material, operator_options_.omega,
                                  operator_options_),
      update_solution_(sph_body, names_.solution, Real(0), names_.search),
      update_residual_(sph_body, names_.residual, Real(1), Real(0), names_.residual, names_.v),
      update_search_(sph_body, names_.search, Real(0), names_.z),
      form_residual_gap_(sph_body, names_.t, Real(1), Real(-1), names_.residual, names_.true_residual),
      dot_r_z_(sph_body, names_.residual, names_.z),
      dot_p_q_(sph_body, names_.search, names_.v),
      norm_squared_r_(sph_body, names_.residual),
      norm_squared_true_r_(sph_body, names_.true_residual),
      norm_squared_gap_(sph_body, names_.t)
{
}

template <class ExecutionPolicy>
inline AphiPCGSolverCK<ExecutionPolicy>::AphiPCGSolverCK(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &operator_options, Real tolerance, UnsignedInt max_iterations,
    bool use_block_jacobi_preconditioner, bool enable_debug_log, Real residual_explosion_factor)
    : AphiPCGSolverCK(sph_body, inner_relation, variable_names, operator_options,
                      AphiPCGSolverOptions{.max_iterations = max_iterations,
                                           .relative_tolerance = tolerance,
                                           .use_block_jacobi_preconditioner = use_block_jacobi_preconditioner,
                                           .enable_debug_log = enable_debug_log,
                                           .residual_explosion_factor = residual_explosion_factor})
{
}

template <class ExecutionPolicy>
inline AphiPCGResult AphiPCGSolverCK<ExecutionPolicy>::finalizeBreakdown(AphiPCGResult &result,
                                                                          AphiPCGBreakdownCode code, Real residual_norm,
                                                                          Real initial_norm)
{
    result.breakdown = true;
    result.breakdown_code = code;
    result.final_residual_norm = residual_norm;
    result.final_relative_residual = residual_norm / initial_norm;
    return result;
}

template <class ExecutionPolicy>
inline void AphiPCGSolverCK<ExecutionPolicy>::maybeRecomputeTrueResidual(UnsignedInt iteration, Real initial_norm,
                                                                        Real &true_relative_residual,
                                                                        Real &recursive_true_gap,
                                                                        Real &true_residual_norm)
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
}

template <class ExecutionPolicy>
inline AphiPCGResult AphiPCGSolverCK<ExecutionPolicy>::solve()
{
    AphiPCGResult result;

    if (solver_options_.use_block_jacobi_preconditioner)
    {
        compute_jacobi_diagonal_.exec();
    }

    apply_solution_to_lhs_.exec();
    compute_recursive_residual_.exec();

    const Real initial_norm = std::sqrt(norm_squared_r_.exec());
    result.initial_residual_norm = initial_norm;
    if (!aphiIsFinite(initial_norm))
    {
        result.breakdown = true;
        result.breakdown_code = AphiPCGBreakdownCode::InitialResidualNonFinite;
        return result;
    }
    if (initial_norm <= solver_options_.absolute_tolerance || initial_norm <= solver_options_.relative_tolerance)
    {
        result.final_residual_norm = initial_norm;
        result.final_relative_residual = Real(0);
        result.converged = true;
        return result;
    }

    if (solver_options_.use_block_jacobi_preconditioner)
    {
        precondition_residual_to_z_.exec();
    }
    else
    {
        copy_r_to_z_.exec();
    }

    copy_z_to_search_.exec();

    Real rho_old = dot_r_z_.exec();
    if (!aphiIsFinite(rho_old) || aphiIsNearZeroAbs(rho_old, solver_options_.coefficient_breakdown_tol))
    {
        return finalizeBreakdown(result, AphiPCGBreakdownCode::RhoNearZero, initial_norm, initial_norm);
    }

    for (UnsignedInt iteration = 0; iteration < solver_options_.max_iterations; ++iteration)
    {
        apply_search_to_v_.exec();

        const Real denom = dot_p_q_.exec();
        if (!aphiIsFinite(denom) || aphiIsNearZeroAbs(denom, solver_options_.coefficient_breakdown_tol))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiPCGBreakdownCode::DenominatorNearZero, residual_norm, initial_norm);
        }
        if (denom <= Real(0))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiPCGBreakdownCode::NonPositiveCurvature, residual_norm, initial_norm);
        }

        const Real alpha = rho_old / denom;
        if (!aphiIsFinite(alpha))
        {
            const Real residual_norm = std::sqrt(norm_squared_r_.exec());
            return finalizeBreakdown(result, AphiPCGBreakdownCode::AlphaNonFinite, residual_norm, initial_norm);
        }

        update_solution_.setAlpha(alpha);
        update_solution_.exec();

        update_residual_.setCoefficients(Real(1), -alpha);
        update_residual_.exec();

        const Real residual_norm = std::sqrt(norm_squared_r_.exec());
        const Real relative_residual = residual_norm / initial_norm;
        result.iteration_count = iteration + 1;
        result.final_residual_norm = residual_norm;
        result.final_relative_residual = relative_residual;

        Real true_relative_residual = relative_residual;
        Real recursive_true_gap = Real(0);
        Real true_residual_norm = result.final_true_residual_norm;
        maybeRecomputeTrueResidual(iteration, initial_norm, true_relative_residual, recursive_true_gap,
                                   true_residual_norm);
        result.final_true_relative_residual = true_relative_residual;
        result.final_recursive_true_gap = recursive_true_gap;
        result.final_true_residual_norm = true_residual_norm;

        if (solver_options_.enable_debug_log)
        {
            std::cout << "[AphiPCG] iter=" << iteration << " alpha=" << alpha << " denom=" << denom
                      << " rho=" << rho_old << " recursive_rel_res=" << relative_residual
                      << " true_rel_res=" << true_relative_residual << " residual_gap=" << recursive_true_gap
                      << std::endl;
        }

        if (!aphiIsFinite(residual_norm))
        {
            return finalizeBreakdown(result, AphiPCGBreakdownCode::ResidualNonFinite, residual_norm, initial_norm);
        }
        if (aphiIsResidualExplosion(relative_residual, solver_options_.residual_explosion_factor))
        {
            return finalizeBreakdown(result, AphiPCGBreakdownCode::ResidualExplosion, residual_norm, initial_norm);
        }

        if (relative_residual < solver_options_.relative_tolerance || residual_norm < solver_options_.absolute_tolerance)
        {
            result.converged = true;
            return result;
        }

        if (solver_options_.use_block_jacobi_preconditioner)
        {
            precondition_residual_to_z_.exec();
        }
        else
        {
            copy_r_to_z_.exec();
        }

        const Real rho_new = dot_r_z_.exec();
        if (!aphiIsFinite(rho_new) || aphiIsNearZeroAbs(rho_new, solver_options_.coefficient_breakdown_tol))
        {
            return finalizeBreakdown(result, AphiPCGBreakdownCode::RhoNearZero, residual_norm, initial_norm);
        }

        if (iteration + 1 < solver_options_.max_iterations)
        {
            const Real beta = rho_new / rho_old;
            if (!aphiIsFinite(beta))
            {
                return finalizeBreakdown(result, AphiPCGBreakdownCode::BetaNonFinite, residual_norm, initial_norm);
            }
            update_search_.setBeta(beta);
            update_search_.exec();
            rho_old = rho_new;
        }
    }

    result.breakdown = true;
    result.breakdown_code = AphiPCGBreakdownCode::MaxIterationsReached;
    return result;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PCG_SOLVER_CK_HPP
