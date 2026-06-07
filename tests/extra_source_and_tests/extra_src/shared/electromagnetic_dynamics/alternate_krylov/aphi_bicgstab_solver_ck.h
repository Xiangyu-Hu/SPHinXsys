#ifndef APHI_BICGSTAB_SOLVER_CK_H
#define APHI_BICGSTAB_SOLVER_CK_H

#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.h"
#include "electromagnetic_dynamics/aphi_block_zero_ck.h"
#include "electromagnetic_dynamics/aphi_block_vector_ops_ck.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.h"
#include "interaction_ck.h"
#include "simple_algorithms_ck.h"

namespace SPH
{
namespace electromagnetics
{

enum class AphiBiCGStabBreakdownCode
{
    None = 0,
    InitialResidualNonFinite,
    RhoOldNearZero,
    RhoNewNearZero,
    RHatDotVNearZero,
    TDotTNearZero,
    OmegaNearZero,
    AlphaNonFinite,
    BetaNonFinite,
    OmegaNonFinite,
    CoefficientNonFinite,
    ResidualNonFinite,
    ResidualExplosion,
    MaxIterationsReached
};

inline const char *AphiBiCGStabBreakdownCodeName(AphiBiCGStabBreakdownCode code)
{
    switch (code)
    {
    case AphiBiCGStabBreakdownCode::None:
        return "None";
    case AphiBiCGStabBreakdownCode::InitialResidualNonFinite:
        return "InitialResidualNonFinite";
    case AphiBiCGStabBreakdownCode::RhoOldNearZero:
        return "RhoOldNearZero";
    case AphiBiCGStabBreakdownCode::RhoNewNearZero:
        return "RhoNewNearZero";
    case AphiBiCGStabBreakdownCode::RHatDotVNearZero:
        return "RHatDotVNearZero";
    case AphiBiCGStabBreakdownCode::TDotTNearZero:
        return "TDotTNearZero";
    case AphiBiCGStabBreakdownCode::OmegaNearZero:
        return "OmegaNearZero";
    case AphiBiCGStabBreakdownCode::AlphaNonFinite:
        return "AlphaNonFinite";
    case AphiBiCGStabBreakdownCode::BetaNonFinite:
        return "BetaNonFinite";
    case AphiBiCGStabBreakdownCode::OmegaNonFinite:
        return "OmegaNonFinite";
    case AphiBiCGStabBreakdownCode::CoefficientNonFinite:
        return "CoefficientNonFinite";
    case AphiBiCGStabBreakdownCode::ResidualNonFinite:
        return "ResidualNonFinite";
    case AphiBiCGStabBreakdownCode::ResidualExplosion:
        return "ResidualExplosion";
    case AphiBiCGStabBreakdownCode::MaxIterationsReached:
        return "MaxIterationsReached";
    default:
        return "Unknown";
    }
}

struct AphiBiCGStabSolverOptions
{
    UnsignedInt max_iterations = 500;
    Real relative_tolerance = Real(1.0e-8);
    Real absolute_tolerance = Real(1.0e-12);

    bool use_block_jacobi_preconditioner = false;
    bool enable_debug_log = false;
    bool recompute_true_residual = false;
    UnsignedInt true_residual_check_interval = 10;

    Real coefficient_breakdown_tol = Real(1.0e-30);
    Real residual_explosion_factor = Real(1.0e12);
    Real denominator_relative_tol = Real(1.0e-14);
};

struct AphiBiCGStabResult
{
    UnsignedInt iteration_count = 0;
    Real initial_residual_norm = 0.0;
    Real final_residual_norm = 0.0;
    Real final_relative_residual = 0.0;
    Real final_true_residual_norm = 0.0;
    Real final_true_relative_residual = 0.0;
    Real final_recursive_true_gap = 0.0;
    bool converged = false;
    bool breakdown = false;
    AphiBiCGStabBreakdownCode breakdown_code = AphiBiCGStabBreakdownCode::None;
};

/**
 * BiCGStab for discrete K(X)=RHS using fused AphiApplyCK.
 * Optional right-preconditioned block-Jacobi:
 *     z = M^{-1} p,   v = A z
 *     y = M^{-1} s,   t = A y
 *     x <- x + alpha z + omega y
 */
template <class ExecutionPolicy>
class AphiBiCGStabSolverCK
{
  public:
    AphiBiCGStabSolverCK(SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
                         const AphiLhsAssemblyOptions &operator_options, const AphiBiCGStabSolverOptions &solver_options);

    AphiBiCGStabSolverCK(SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
                         const AphiLhsAssemblyOptions &operator_options, Real tolerance, UnsignedInt max_iterations,
                         bool use_block_jacobi_preconditioner = false, bool enable_debug_log = false,
                         Real residual_explosion_factor = Real(1.0e12));

    AphiBiCGStabResult solve();

  protected:
    AphiBiCGStabResult finalizeBreakdown(AphiBiCGStabResult &result, AphiBiCGStabBreakdownCode code, Real residual_norm,
                                         Real initial_norm);

    void maybeRecomputeTrueResidual(UnsignedInt iteration, Real initial_norm, Real recursive_norm,
                                    Real &true_relative_residual, Real &recursive_true_gap, Real &true_residual_norm);

    SPHBody &body_;
    AphiVariableNames names_;
    AphiLhsAssemblyOptions operator_options_;
    AphiBiCGStabSolverOptions solver_options_;

    InteractionDynamicsCK<ExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_solution_to_lhs_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_search_to_v_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_z_to_v_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_y_to_v_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_s_to_t_;
    StateDynamics<ExecutionPolicy, AphiComputeResidualCK> compute_recursive_residual_;
    StateDynamics<ExecutionPolicy, AphiComputeBlockResidualCK> compute_true_residual_;
    StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_residual_to_r_hat_;
    StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_residual_to_search_;
    StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_v_to_v_old_;
    StateDynamics<ExecutionPolicy, AphiApplyBlockJacobiInverseCK> precondition_search_to_z_;
    StateDynamics<ExecutionPolicy, AphiApplyBlockJacobiInverseCK> precondition_s_to_y_;
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_s_;
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_residual_gap_;
    StateDynamics<ExecutionPolicy, AphiBlockBiCGStabUpdateSolutionCK> update_solution_;
    StateDynamics<ExecutionPolicy, AphiBlockBiCGStabUpdateSolutionCK> update_solution_preconditioned_;
    StateDynamics<ExecutionPolicy, AphiBlockBiCGStabUpdateResidualCK> update_residual_;
    StateDynamics<ExecutionPolicy, AphiBlockBiCGStabUpdateResidualCK> update_residual_preconditioned_;
    StateDynamics<ExecutionPolicy, AphiBlockBiCGStabUpdateSearchCK> update_search_;
    StateDynamics<ExecutionPolicy, AphiBlockBiCGStabUpdateSearchCK> update_search_preconditioned_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_r_hat_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_r_hat_v_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_v_s_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_t_s_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_v_v_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_t_t_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_true_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_gap_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_BICGSTAB_SOLVER_CK_H
