#ifndef APHI_ASSEMBLE_LHS_DEBUG_CK_H
#define APHI_ASSEMBLE_LHS_DEBUG_CK_H

#include "interaction_ck.h"
#include "electromagnetic_dynamics/aphi_block_zero_ck.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_div_sigma_a_coupling_ck.h"
#include "electromagnetic_dynamics/aphi_grad_phi_coupling_ck.h"
#include "electromagnetic_dynamics/aphi_gradient_divergence_debug_ck.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"
#include "electromagnetic_dynamics/aphi_reaction_ck.h"
#include "general_gradient.h"
#include "kernel_correction_ck.h"

namespace SPH
{
namespace electromagnetics
{

/**
 * Stage 3B debug LHS assembler: K(X) via composable terms.
 * OPERATOR default uses PairwiseUncorrected grad/div modes.
 */
template <class ExecutionPolicy>
class AphiAssembleLhsDebugDynamicsBundle
{
  public:
    AphiAssembleLhsDebugDynamicsBundle(SPHBody &sph_body, Inner<> &inner_relation,
                                       const AphiVariableNames &variable_names,
                                       const AphiLhsAssemblyOptions &options);

    void exec();

  protected:
    AphiLhsAssemblyOptions options_;
    AphiVariableNames names_;

    StateDynamics<ExecutionPolicy, AphiZeroBlockCK> zero_lhs_;
    StateDynamics<ExecutionPolicy, AphiReactionCK> reaction_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>>> laplace_phi_real_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>>> laplace_phi_imag_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Vecd>>> laplace_a_real_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Vecd>>> laplace_a_imag_;
    InteractionDynamicsCK<ExecutionPolicy, AphiGradPhiCouplingCK<Inner<>>> grad_phi_pairwise_;
    InteractionDynamicsCK<ExecutionPolicy, AphiDivSigmaACouplingCK<Inner<>>> div_sigma_a_pairwise_;

    InteractionDynamicsCK<ExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Real>>> phi_real_gradient_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Real>>> phi_imag_gradient_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient_;
    StateDynamics<ExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_real_diagnostic_;
    StateDynamics<ExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_imag_diagnostic_;
    StateDynamics<ExecutionPolicy, AphiGradPhiDiagnosticCouplingCK> grad_phi_diagnostic_;
    StateDynamics<ExecutionPolicy, AphiDivSigmaAConstSigmaDiagnosticCouplingCK> div_sigma_a_diagnostic_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_ASSEMBLE_LHS_DEBUG_CK_H
