#ifndef APHI_MATRIX_FREE_SOLVE_CK_H
#define APHI_MATRIX_FREE_SOLVE_CK_H

#include "electromagnetic_dynamics/alternate_krylov/aphi_bicgstab_solver_ck.h"
#include "electromagnetic_dynamics/aphi_gmres_solver_ck.h"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.h"
#include "electromagnetic_dynamics/alternate_krylov/aphi_pcg_solver_ck.h"

namespace SPH
{
namespace electromagnetics
{

enum class AphiKrylovSolverKind
{
    /** Production default: restarted right-preconditioned GMRES(m). */
    GMRES,
    /** Diagnostic only: BiCGStab with optional block-Jacobi PC. */
    BiCGStabDiagnostic,
    /** Diagnostic only: scalar block PCG (SPD assumption). */
    PCGScalarDiagnostic
};

inline const char *AphiKrylovSolverKindName(AphiKrylovSolverKind kind)
{
    switch (kind)
    {
    case AphiKrylovSolverKind::GMRES:
        return "GMRES";
    case AphiKrylovSolverKind::BiCGStabDiagnostic:
        return "BiCGStabDiagnostic";
    case AphiKrylovSolverKind::PCGScalarDiagnostic:
        return "PCGScalarDiagnostic";
    default:
        return "Unknown";
    }
}

struct AphiMatrixFreeSolverOptions
{
    AphiKrylovSolverKind solver_kind = AphiKrylovSolverKind::GMRES;
    AphiGMRESSolverOptions gmres{};
    AphiBiCGStabSolverOptions bicgstab{};
    AphiPCGSolverOptions pcg{};
};

struct AphiMatrixFreeSolverResult
{
    AphiKrylovSolverKind solver_kind = AphiKrylovSolverKind::GMRES;
    UnsignedInt outer_iteration_count = 0;
    UnsignedInt arnoldi_step_count = 0;
    UnsignedInt iteration_count = 0;
    Real initial_residual_norm = 0.0;
    Real final_relative_residual = 0.0;
    Real final_true_relative_residual = 0.0;
    Real final_true_residual_norm = 0.0;
    Real final_recursive_true_gap = 0.0;
    Real final_res_a_total_norm = 0.0;
    Real final_res_phi_total_norm = 0.0;
    Real final_res_a_fraction = 0.0;
    Real final_res_phi_fraction = 0.0;
    bool monotonic_outer_residual = true;
    bool converged = false;
    bool breakdown = false;
    const char *breakdown_code_name = "None";
};

/**
 * High-level matrix-free Krylov dispatch. Production path defaults to GMRES.
 * GMRES workspace (GMRESV/GMRESZ blocks) is registered in the constructor.
 */
template <class ExecutionPolicy>
class AphiMatrixFreeSolveCK
{
  public:
    AphiMatrixFreeSolveCK(SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
                          const AphiLhsAssemblyOptions &operator_options,
                          const AphiMatrixFreeSolverOptions &solver_options);

    AphiMatrixFreeSolverResult solve();

  protected:
    AphiMatrixFreeSolverResult fromGMRESResult(const AphiGMRESResult &result) const;
    AphiMatrixFreeSolverResult fromBiCGStabResult(const AphiBiCGStabResult &result) const;
    AphiMatrixFreeSolverResult fromPCGResult(const AphiPCGResult &result) const;

    SPHBody &body_;
    Inner<> &inner_;
    AphiVariableNames names_;
    AphiLhsAssemblyOptions operator_options_;
    AphiMatrixFreeSolverOptions solver_options_;
    AphiGMRESWorkspaceNames gmres_workspace_;
    RegisterAphiGMRESWorkspaceCK register_gmres_workspace_;
};

/** Matrix-free Krylov dispatch for Inner+Contact multi-body relations (GMRES only). */
template <class ExecutionPolicy>
class AphiMatrixFreeContactSolveCK
{
  public:
    AphiMatrixFreeContactSolveCK(SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation,
                                 const AphiVariableNames &variable_names, const AphiLhsAssemblyOptions &operator_options,
                                 const AphiMatrixFreeSolverOptions &solver_options);

    AphiMatrixFreeSolverResult solve();

  protected:
    AphiMatrixFreeSolverResult fromGMRESResult(const AphiGMRESResult &result) const;

    SPHBody &body_;
    Inner<> &inner_;
    Contact<> &contact_;
    AphiVariableNames names_;
    AphiLhsAssemblyOptions operator_options_;
    AphiMatrixFreeSolverOptions solver_options_;
    AphiGMRESWorkspaceNames gmres_workspace_;
    RegisterAphiGMRESWorkspaceCK register_gmres_workspace_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_MATRIX_FREE_SOLVE_CK_H
