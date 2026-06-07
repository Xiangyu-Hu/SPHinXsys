#ifndef APHI_GMRES_SOLVER_CK_H
#define APHI_GMRES_SOLVER_CK_H

#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.h"
#include "electromagnetic_dynamics/aphi_block_vector_ops_ck.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.h"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
namespace electromagnetics
{

enum class AphiGMRESBreakdownCode
{
    None = 0,
    InitialResidualNonFinite,
    HappyBreakdown,
    ArnoldiNormNonFinite,
    LeastSquaresFailed,
    ResidualNonFinite,
    ResidualExplosion,
    MaxOuterIterationsReached
};

inline const char *AphiGMRESBreakdownCodeName(AphiGMRESBreakdownCode code)
{
    switch (code)
    {
    case AphiGMRESBreakdownCode::None:
        return "None";
    case AphiGMRESBreakdownCode::InitialResidualNonFinite:
        return "InitialResidualNonFinite";
    case AphiGMRESBreakdownCode::HappyBreakdown:
        return "HappyBreakdown";
    case AphiGMRESBreakdownCode::ArnoldiNormNonFinite:
        return "ArnoldiNormNonFinite";
    case AphiGMRESBreakdownCode::LeastSquaresFailed:
        return "LeastSquaresFailed";
    case AphiGMRESBreakdownCode::ResidualNonFinite:
        return "ResidualNonFinite";
    case AphiGMRESBreakdownCode::ResidualExplosion:
        return "ResidualExplosion";
    case AphiGMRESBreakdownCode::MaxOuterIterationsReached:
        return "MaxOuterIterationsReached";
    default:
        return "Unknown";
    }
}

enum class AphiMultiBodyInnerProductKind
{
    /** Sum vol-weighted contributions from all bodies (particle-count / volume dominated). */
    VolWeightedGlobalSum,
    /** Average per-body vol-weighted inner products (balanced across bodies). */
    EqualPerBodyAverage,
};

struct AphiGMRESSolverOptions
{
    UnsignedInt restart_dimension = 30;
    UnsignedInt max_outer_iterations = 20;
    Real relative_tolerance = Real(1.0e-8);
    Real absolute_tolerance = Real(1.0e-12);
    AphiMultiBodyInnerProductKind multibody_inner_product = AphiMultiBodyInnerProductKind::EqualPerBodyAverage;

    bool use_block_jacobi_preconditioner = true;
    AphiBlockJacobiPreconditionerKind block_jacobi_kind = AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8;
    bool enable_debug_log = false;
    bool recompute_true_residual = true;
    UnsignedInt true_residual_check_interval = 5;

    Real happy_breakdown_tol = Real(1.0e-14);
    Real residual_explosion_factor = Real(1.0e12);
};

struct AphiGMRESResult
{
    UnsignedInt outer_iteration_count = 0;
    UnsignedInt arnoldi_step_count = 0;
    Real initial_residual_norm = 0.0;
    Real final_residual_norm = 0.0;
    Real final_relative_residual = 0.0;
    Real final_true_residual_norm = 0.0;
    Real final_true_relative_residual = 0.0;
    Real final_recursive_true_gap = 0.0;
    Real final_res_a_total_norm = 0.0;
    Real final_res_phi_total_norm = 0.0;
    Real final_res_a_fraction = 0.0;
    Real final_res_phi_fraction = 0.0;
    bool monotonic_outer_residual = true;
    bool converged = false;
    bool breakdown = false;
    AphiGMRESBreakdownCode breakdown_code = AphiGMRESBreakdownCode::None;
};

/**
 * Restarted right-preconditioned GMRES for discrete K(X)=RHS.
 * Production default Krylov backend for A-phi matrix-free solves.
 * Workspace blocks: GMRESV0..GMRESV{m}, GMRESZ0..GMRESZ{m-1}.
 * Algorithm: z_j = M^{-1} v_j, w = A z_j, x += sum y_j z_j.
 */
template <class ExecutionPolicy>
class AphiGMRESSolverCK
{
  public:
    AphiGMRESSolverCK(SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
                      const AphiGMRESWorkspaceNames &workspace_names, const AphiLhsAssemblyOptions &operator_options,
                      const AphiGMRESSolverOptions &solver_options);

    AphiGMRESResult solve();

  protected:
    SPHBody &body_;
    Inner<> &inner_;
    AphiVariableNames names_;
    AphiGMRESWorkspaceNames workspace_;
    AphiLhsAssemblyOptions operator_options_;
    AphiGMRESSolverOptions solver_options_;

    InteractionDynamicsCK<ExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_solution_to_lhs_;
    StateDynamics<ExecutionPolicy, AphiComputeResidualCK> compute_recursive_residual_;
    StateDynamics<ExecutionPolicy, AphiComputeBlockResidualCK> compute_true_residual_;
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_residual_gap_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_true_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_gap_;
};

/** Restarted right-preconditioned GMRES with Inner+Contact fused apply and block-Jacobi diagonal. */
template <class ExecutionPolicy>
class AphiGMRESContactSolverCK
{
  public:
    AphiGMRESContactSolverCK(SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation,
                             const AphiVariableNames &variable_names, const AphiGMRESWorkspaceNames &workspace_names,
                             const AphiLhsAssemblyOptions &operator_options, const AphiGMRESSolverOptions &solver_options);

    AphiGMRESResult solve();

  protected:
    SPHBody &body_;
    Inner<> &inner_;
    Contact<> &contact_;
    AphiVariableNames names_;
    AphiGMRESWorkspaceNames workspace_;
    AphiLhsAssemblyOptions operator_options_;
    AphiGMRESSolverOptions solver_options_;

    AphiComputeBlockJacobiContactDynamicsBundle<ExecutionPolicy> compute_jacobi_diagonal_;
    AphiApplyContactDynamicsBundle<ExecutionPolicy> apply_solution_to_lhs_;
    StateDynamics<ExecutionPolicy, AphiComputeResidualCK> compute_recursive_residual_;
    StateDynamics<ExecutionPolicy, AphiComputeBlockResidualCK> compute_true_residual_;
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_residual_gap_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_true_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_gap_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GMRES_SOLVER_CK_H
