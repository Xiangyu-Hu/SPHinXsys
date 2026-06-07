#ifndef APHI_GMRES_TEST_HELPERS_H
#define APHI_GMRES_TEST_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.hpp"
#include "electromagnetic_dynamics/aphi_matrix_free_solve_ck.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline AphiGMRESSolverOptions defaultGMRESConvergenceOptions(Real tolerance, UnsignedInt restart_dimension = 50,
                                                             UnsignedInt max_outer_iterations = 20)
{
    AphiGMRESSolverOptions options;
    options.restart_dimension = restart_dimension;
    options.max_outer_iterations = max_outer_iterations;
    options.relative_tolerance = tolerance;
    options.use_block_jacobi_preconditioner = true;
    options.block_jacobi_kind = AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8;
    options.recompute_true_residual = true;
    options.true_residual_check_interval = 5;
    options.enable_debug_log = false;
    return options;
}

inline AphiMatrixFreeSolverOptions defaultMatrixFreeGMRESOptions(Real tolerance, UnsignedInt restart_dimension = 50,
                                                                 UnsignedInt max_outer_iterations = 20)
{
    AphiMatrixFreeSolverOptions options;
    options.solver_kind = AphiKrylovSolverKind::GMRES;
    options.gmres = defaultGMRESConvergenceOptions(tolerance, restart_dimension, max_outer_iterations);
    return options;
}

/** Coupled multi-body Contact GMRES defaults: equal-per-body Krylov inner products. */
inline AphiMatrixFreeSolverOptions defaultCoupledContactGMRESOptions(Real tolerance,
                                                                     UnsignedInt restart_dimension = 50,
                                                                     UnsignedInt max_outer_iterations = 100)
{
    AphiMatrixFreeSolverOptions options = defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);
    options.gmres.multibody_inner_product = AphiMultiBodyInnerProductKind::EqualPerBodyAverage;
    return options;
}

inline bool gmresConvergencePassed(const AphiGMRESResult &result, Real tolerance, Real true_tolerance_factor = 10.0,
                                   Real max_recursive_true_gap = 1.0e-4)
{
    return result.converged && !result.breakdown && result.final_relative_residual < tolerance &&
           result.final_true_relative_residual < true_tolerance_factor * tolerance &&
           result.final_recursive_true_gap < max_recursive_true_gap;
}

inline bool gmresConvergencePassed(const AphiMatrixFreeSolverResult &result, Real tolerance,
                                   Real true_tolerance_factor = 10.0, Real max_recursive_true_gap = 1.0e-4)
{
    return result.solver_kind == AphiKrylovSolverKind::GMRES && result.converged && !result.breakdown &&
           result.final_relative_residual < tolerance &&
           result.final_true_relative_residual < true_tolerance_factor * tolerance &&
           result.final_recursive_true_gap < max_recursive_true_gap;
}

inline AphiMatrixFreeSolverResult toMatrixFreeSolverResult(const AphiGMRESResult &result)
{
    AphiMatrixFreeSolverResult unified;
    unified.solver_kind = AphiKrylovSolverKind::GMRES;
    unified.outer_iteration_count = result.outer_iteration_count;
    unified.arnoldi_step_count = result.arnoldi_step_count;
    unified.iteration_count = result.arnoldi_step_count;
    unified.initial_residual_norm = result.initial_residual_norm;
    unified.final_relative_residual = result.final_relative_residual;
    unified.final_true_relative_residual = result.final_true_relative_residual;
    unified.final_true_residual_norm = result.final_true_residual_norm;
    unified.final_recursive_true_gap = result.final_recursive_true_gap;
    unified.monotonic_outer_residual = result.monotonic_outer_residual;
    unified.converged = result.converged;
    unified.breakdown = result.breakdown;
    unified.breakdown_code_name = AphiGMRESBreakdownCodeName(result.breakdown_code);
    return unified;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_TEST_HELPERS_H
