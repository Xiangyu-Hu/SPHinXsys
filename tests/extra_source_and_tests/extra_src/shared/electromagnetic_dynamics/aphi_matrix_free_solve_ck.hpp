#ifndef APHI_MATRIX_FREE_SOLVE_CK_HPP
#define APHI_MATRIX_FREE_SOLVE_CK_HPP

#include "electromagnetic_dynamics/aphi_matrix_free_solve_ck.h"
#include "electromagnetic_dynamics/alternate_krylov/aphi_bicgstab_solver_ck.hpp"
#include "electromagnetic_dynamics/alternate_krylov/aphi_pcg_solver_ck.hpp"
#include "electromagnetic_dynamics/aphi_gmres_solver_ck.hpp"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.hpp"

namespace SPH
{
namespace electromagnetics
{

template <class ExecutionPolicy>
inline AphiMatrixFreeSolveCK<ExecutionPolicy>::AphiMatrixFreeSolveCK(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &operator_options, const AphiMatrixFreeSolverOptions &solver_options)
    : body_(sph_body), inner_(inner_relation), names_(variable_names), operator_options_(operator_options),
      solver_options_(solver_options),
      gmres_workspace_(buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension)),
      register_gmres_workspace_(sph_body, solver_options.gmres.restart_dimension)
{
    (void)register_gmres_workspace_;
}

template <class ExecutionPolicy>
inline AphiMatrixFreeSolverResult
AphiMatrixFreeSolveCK<ExecutionPolicy>::fromGMRESResult(const AphiGMRESResult &result) const
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
    unified.final_res_a_total_norm = result.final_res_a_total_norm;
    unified.final_res_phi_total_norm = result.final_res_phi_total_norm;
    unified.final_res_a_fraction = result.final_res_a_fraction;
    unified.final_res_phi_fraction = result.final_res_phi_fraction;
    unified.monotonic_outer_residual = result.monotonic_outer_residual;
    unified.converged = result.converged;
    unified.breakdown = result.breakdown;
    unified.breakdown_code_name = AphiGMRESBreakdownCodeName(result.breakdown_code);
    return unified;
}

template <class ExecutionPolicy>
inline AphiMatrixFreeSolverResult
AphiMatrixFreeSolveCK<ExecutionPolicy>::fromBiCGStabResult(const AphiBiCGStabResult &result) const
{
    AphiMatrixFreeSolverResult unified;
    unified.solver_kind = AphiKrylovSolverKind::BiCGStabDiagnostic;
    unified.iteration_count = result.iteration_count;
    unified.initial_residual_norm = result.initial_residual_norm;
    unified.final_relative_residual = result.final_relative_residual;
    unified.final_true_relative_residual = result.final_true_relative_residual;
    unified.final_true_residual_norm = result.final_true_residual_norm;
    unified.final_recursive_true_gap = result.final_recursive_true_gap;
    unified.converged = result.converged;
    unified.breakdown = result.breakdown;
    unified.breakdown_code_name = AphiBiCGStabBreakdownCodeName(result.breakdown_code);
    return unified;
}

template <class ExecutionPolicy>
inline AphiMatrixFreeSolverResult
AphiMatrixFreeSolveCK<ExecutionPolicy>::fromPCGResult(const AphiPCGResult &result) const
{
    AphiMatrixFreeSolverResult unified;
    unified.solver_kind = AphiKrylovSolverKind::PCGScalarDiagnostic;
    unified.iteration_count = result.iteration_count;
    unified.initial_residual_norm = result.initial_residual_norm;
    unified.final_relative_residual = result.final_relative_residual;
    unified.final_true_relative_residual = result.final_true_relative_residual;
    unified.final_true_residual_norm = result.final_true_residual_norm;
    unified.final_recursive_true_gap = result.final_recursive_true_gap;
    unified.converged = result.converged;
    unified.breakdown = result.breakdown;
    unified.breakdown_code_name = AphiPCGBreakdownCodeName(result.breakdown_code);
    return unified;
}

template <class ExecutionPolicy>
inline AphiMatrixFreeSolverResult AphiMatrixFreeSolveCK<ExecutionPolicy>::solve()
{
    switch (solver_options_.solver_kind)
    {
    case AphiKrylovSolverKind::BiCGStabDiagnostic:
    {
        AphiBiCGStabSolverCK<ExecutionPolicy> solver(body_, inner_, names_, operator_options_,
                                                     solver_options_.bicgstab);
        return fromBiCGStabResult(solver.solve());
    }
    case AphiKrylovSolverKind::PCGScalarDiagnostic:
    {
        AphiPCGSolverCK<ExecutionPolicy> solver(body_, inner_, names_, operator_options_, solver_options_.pcg);
        return fromPCGResult(solver.solve());
    }
    case AphiKrylovSolverKind::GMRES:
    default:
    {
        AphiGMRESSolverCK<ExecutionPolicy> solver(body_, inner_, names_, gmres_workspace_, operator_options_,
                                                  solver_options_.gmres);
        return fromGMRESResult(solver.solve());
    }
    }
}

template <class ExecutionPolicy>
inline AphiMatrixFreeContactSolveCK<ExecutionPolicy>::AphiMatrixFreeContactSolveCK(
    SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &operator_options, const AphiMatrixFreeSolverOptions &solver_options)
    : body_(sph_body), inner_(inner_relation), contact_(contact_relation), names_(variable_names),
      operator_options_(operator_options), solver_options_(solver_options),
      gmres_workspace_(buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension)),
      register_gmres_workspace_(sph_body, solver_options.gmres.restart_dimension)
{
    (void)register_gmres_workspace_;
}

template <class ExecutionPolicy>
inline AphiMatrixFreeSolverResult
AphiMatrixFreeContactSolveCK<ExecutionPolicy>::fromGMRESResult(const AphiGMRESResult &result) const
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
    unified.final_res_a_total_norm = result.final_res_a_total_norm;
    unified.final_res_phi_total_norm = result.final_res_phi_total_norm;
    unified.final_res_a_fraction = result.final_res_a_fraction;
    unified.final_res_phi_fraction = result.final_res_phi_fraction;
    unified.monotonic_outer_residual = result.monotonic_outer_residual;
    unified.converged = result.converged;
    unified.breakdown = result.breakdown;
    unified.breakdown_code_name = AphiGMRESBreakdownCodeName(result.breakdown_code);
    return unified;
}

template <class ExecutionPolicy>
inline AphiMatrixFreeSolverResult AphiMatrixFreeContactSolveCK<ExecutionPolicy>::solve()
{
    AphiGMRESContactSolverCK<ExecutionPolicy> solver(body_, inner_, contact_, names_, gmres_workspace_,
                                                   operator_options_, solver_options_.gmres);
    return fromGMRESResult(solver.solve());
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_MATRIX_FREE_SOLVE_CK_HPP
