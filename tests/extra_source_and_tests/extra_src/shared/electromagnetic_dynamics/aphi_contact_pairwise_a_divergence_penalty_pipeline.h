#ifndef APHI_CONTACT_PAIRWISE_A_DIVERGENCE_PENALTY_PIPELINE_H
#define APHI_CONTACT_PAIRWISE_A_DIVERGENCE_PENALTY_PIPELINE_H

#include "electromagnetic_dynamics/aphi_a_divergence_penalty_ck.hpp"
#include "electromagnetic_dynamics/aphi_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_pairwise_div_a_ck.hpp"

namespace SPH
{
namespace electromagnetics
{

/** Inner+Contact pairwise divA -> grad(divA) -> LhsA -= lambda*grad(divA). */
template <class ExecutionPolicy>
class AphiContactPairwiseADivergencePenaltyPipelineBundle
{
  public:
    AphiContactPairwiseADivergencePenaltyPipelineBundle(
        SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation, const AphiBlockNames &input_block,
        const AphiBlockNames &output_block, Real a_divergence_penalty,
        const AphiADivergencePenaltyScratchNames &scratch = aphiDefaultADivergencePenaltyScratchNames())
        : a_divergence_penalty_(a_divergence_penalty),
          div_a_real_(DynamicsArgs(inner_relation, input_block.a_real, scratch.div_a_real),
                      DynamicsArgs(contact_relation, input_block.a_real, scratch.div_a_real)),
          div_a_imag_(DynamicsArgs(inner_relation, input_block.a_imag, scratch.div_a_imag),
                      DynamicsArgs(contact_relation, input_block.a_imag, scratch.div_a_imag)),
          div_a_real_gradient_(DynamicsArgs(inner_relation, scratch.div_a_real, scratch.grad_div_a_real),
                               DynamicsArgs(contact_relation, scratch.div_a_real, scratch.grad_div_a_real)),
          div_a_imag_gradient_(DynamicsArgs(inner_relation, scratch.div_a_imag, scratch.grad_div_a_imag),
                               DynamicsArgs(contact_relation, scratch.div_a_imag, scratch.grad_div_a_imag)),
          grad_div_penalty_(sph_body, scratch.grad_div_a_real, scratch.grad_div_a_imag, output_block,
                            a_divergence_penalty)
    {
    }

    void exec()
    {
        div_a_real_.exec();
        div_a_imag_.exec();
        div_a_real_gradient_.exec();
        div_a_imag_gradient_.exec();
        grad_div_penalty_.exec();
    }

  protected:
    Real a_divergence_penalty_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseVectorDivergenceCK<Inner<>, Contact<>>> div_a_real_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseVectorDivergenceCK<Inner<>, Contact<>>> div_a_imag_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseScalarGradientCK<Inner<>, Contact<>>> div_a_real_gradient_;
    InteractionDynamicsCK<ExecutionPolicy, AphiPairwiseScalarGradientCK<Inner<>, Contact<>>> div_a_imag_gradient_;
    StateDynamics<ExecutionPolicy, AphiGradDivAPenaltyCK> grad_div_penalty_;
};

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_PAIRWISE_A_DIVERGENCE_PENALTY_PIPELINE_H
