#ifndef APHI_A_DIVERGENCE_PENALTY_PIPELINE_H
#define APHI_A_DIVERGENCE_PENALTY_PIPELINE_H

#include "electromagnetic_dynamics/aphi_a_divergence_penalty_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_gradient_divergence_debug_ck.hpp"
#include "general_gradient.h"
#include "kernel_correction_ck.h"

namespace SPH
{
namespace electromagnetics
{

/** B-corrected grad(div A) penalty contribution: output A += lambda * grad(div A) from input A block. */
template <class ExecutionPolicy>
class AphiADivergencePenaltyPipelineBundle
{
  public:
    AphiADivergencePenaltyPipelineBundle(SPHBody &sph_body, Inner<> &inner_relation, const AphiBlockNames &input_block,
                                         const AphiBlockNames &output_block, Real a_divergence_penalty,
                                         const AphiADivergencePenaltyScratchNames &scratch =
                                             aphiDefaultADivergencePenaltyScratchNames())
        : a_divergence_penalty_(a_divergence_penalty),
          linear_correction_matrix_(DynamicsArgs(inner_relation, 0.0)),
          a_real_gradient_(DynamicsArgs(inner_relation, input_block.a_real)),
          a_imag_gradient_(DynamicsArgs(inner_relation, input_block.a_imag)),
          div_a_real_(sph_body, input_block.a_real + "Gradient", scratch.div_a_real),
          div_a_imag_(sph_body, input_block.a_imag + "Gradient", scratch.div_a_imag),
          div_a_real_gradient_(DynamicsArgs(inner_relation, scratch.div_a_real)),
          div_a_imag_gradient_(DynamicsArgs(inner_relation, scratch.div_a_imag)),
          grad_div_penalty_(sph_body, scratch.grad_div_a_real, scratch.grad_div_a_imag, output_block,
                            a_divergence_penalty)
    {
    }

    void exec()
    {
        linear_correction_matrix_.exec();
        a_real_gradient_.exec();
        a_imag_gradient_.exec();
        div_a_real_.exec();
        div_a_imag_.exec();
        div_a_real_gradient_.exec();
        div_a_imag_gradient_.exec();
        grad_div_penalty_.exec();
    }

  protected:
    Real a_divergence_penalty_;
    InteractionDynamicsCK<ExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient_;
    StateDynamics<ExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_real_;
    StateDynamics<ExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_imag_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Real>>> div_a_real_gradient_;
    InteractionDynamicsCK<ExecutionPolicy, LinearGradient<Inner<Real>>> div_a_imag_gradient_;
    StateDynamics<ExecutionPolicy, AphiGradDivAPenaltyCK> grad_div_penalty_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_A_DIVERGENCE_PENALTY_PIPELINE_H
