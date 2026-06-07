#ifndef APHI_PCG_SOLVER_CK_H
#define APHI_PCG_SOLVER_CK_H

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

enum class AphiPCGBreakdownCode
{
    None = 0,
    InitialResidualNonFinite,
    RhoNearZero,
    DenominatorNearZero,
    NonPositiveCurvature,
    AlphaNonFinite,
    BetaNonFinite,
    CoefficientNonFinite,
    ResidualNonFinite,
    ResidualExplosion,
    MaxIterationsReached
};

inline const char *AphiPCGBreakdownCodeName(AphiPCGBreakdownCode code)
{
    switch (code)
    {
    case AphiPCGBreakdownCode::None:
        return "None";
    case AphiPCGBreakdownCode::InitialResidualNonFinite:
        return "InitialResidualNonFinite";
    case AphiPCGBreakdownCode::RhoNearZero:
        return "RhoNearZero";
    case AphiPCGBreakdownCode::DenominatorNearZero:
        return "DenominatorNearZero";
    case AphiPCGBreakdownCode::NonPositiveCurvature:
        return "NonPositiveCurvature";
    case AphiPCGBreakdownCode::AlphaNonFinite:
        return "AlphaNonFinite";
    case AphiPCGBreakdownCode::BetaNonFinite:
        return "BetaNonFinite";
    case AphiPCGBreakdownCode::CoefficientNonFinite:
        return "CoefficientNonFinite";
    case AphiPCGBreakdownCode::ResidualNonFinite:
        return "ResidualNonFinite";
    case AphiPCGBreakdownCode::ResidualExplosion:
        return "ResidualExplosion";
    case AphiPCGBreakdownCode::MaxIterationsReached:
        return "MaxIterationsReached";
    default:
        return "Unknown";
    }
}

struct AphiPCGSolverOptions
{
    UnsignedInt max_iterations = 500;
    Real relative_tolerance = Real(1.0e-8);
    Real absolute_tolerance = Real(1.0e-12);

    bool use_block_jacobi_preconditioner = true;
    bool enable_debug_log = false;
    bool recompute_true_residual = false;
    UnsignedInt true_residual_check_interval = 10;

    Real coefficient_breakdown_tol = Real(1.0e-30);
    Real residual_explosion_factor = Real(1.0e12);
};

struct AphiPCGResult
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
    AphiPCGBreakdownCode breakdown_code = AphiPCGBreakdownCode::None;
};

/** Preconditioned CG for discrete K(X)=RHS using fused AphiApplyCK. */
template <class ExecutionPolicy>
class AphiPCGSolverCK
{
  public:
    AphiPCGSolverCK(SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
                    const AphiLhsAssemblyOptions &operator_options, const AphiPCGSolverOptions &solver_options);

    AphiPCGSolverCK(SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
                    const AphiLhsAssemblyOptions &operator_options, Real tolerance, UnsignedInt max_iterations,
                    bool use_block_jacobi_preconditioner = true, bool enable_debug_log = false,
                    Real residual_explosion_factor = Real(1.0e12));

    AphiPCGResult solve();

  protected:
    AphiPCGResult finalizeBreakdown(AphiPCGResult &result, AphiPCGBreakdownCode code, Real residual_norm,
                                    Real initial_norm);

    void maybeRecomputeTrueResidual(UnsignedInt iteration, Real initial_norm, Real &true_relative_residual,
                                    Real &recursive_true_gap, Real &true_residual_norm);

    SPHBody &body_;
    AphiVariableNames names_;
    AphiLhsAssemblyOptions operator_options_;
    AphiPCGSolverOptions solver_options_;

    InteractionDynamicsCK<ExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_solution_to_lhs_;
    AphiApplyDynamicsBundle<ExecutionPolicy> apply_search_to_v_;
    StateDynamics<ExecutionPolicy, AphiComputeResidualCK> compute_recursive_residual_;
    StateDynamics<ExecutionPolicy, AphiComputeBlockResidualCK> compute_true_residual_;
    StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_z_to_search_;
    StateDynamics<ExecutionPolicy, AphiCopyBlockCK> copy_r_to_z_;
    StateDynamics<ExecutionPolicy, AphiApplyBlockJacobiInverseCK> precondition_residual_to_z_;
    StateDynamics<ExecutionPolicy, AphiBlockAXPYCK> update_solution_;
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> update_residual_;
    StateDynamics<ExecutionPolicy, AphiBlockCGUpdateSearchCK> update_search_;
    StateDynamics<ExecutionPolicy, AphiBlockLinearCombinationCK> form_residual_gap_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_r_z_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockDotProductCK> dot_p_q_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_true_r_;
    ReduceDynamicsCK<ExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_gap_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PCG_SOLVER_CK_H
