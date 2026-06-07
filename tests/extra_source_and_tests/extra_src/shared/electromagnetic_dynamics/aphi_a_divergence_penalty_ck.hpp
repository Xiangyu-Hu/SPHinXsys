#ifndef APHI_A_DIVERGENCE_PENALTY_CK_HPP
#define APHI_A_DIVERGENCE_PENALTY_CK_HPP

#include "electromagnetic_dynamics/aphi_a_divergence_penalty_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiGradDivAPenaltyCK::AphiGradDivAPenaltyCK(SPHBody &sph_body, const std::string &grad_div_a_real_name,
                                                    const std::string &grad_div_a_imag_name,
                                                    const AphiBlockNames &output_block, Real a_divergence_penalty)
    : LocalDynamics(sph_body), a_divergence_penalty_(a_divergence_penalty),
      dv_grad_div_a_real_(particles_->template getVariableByName<Vecd>(grad_div_a_real_name)),
      dv_grad_div_a_imag_(particles_->template getVariableByName<Vecd>(grad_div_a_imag_name)),
      dv_out_a_real_(particles_->template getVariableByName<Vecd>(output_block.a_real)),
      dv_out_a_imag_(particles_->template getVariableByName<Vecd>(output_block.a_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiGradDivAPenaltyCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : a_divergence_penalty_(encloser.a_divergence_penalty_),
      grad_div_a_real_(encloser.dv_grad_div_a_real_->DelegatedData(ex_policy)),
      grad_div_a_imag_(encloser.dv_grad_div_a_imag_->DelegatedData(ex_policy)),
      out_a_real_(encloser.dv_out_a_real_->DelegatedData(ex_policy)),
      out_a_imag_(encloser.dv_out_a_imag_->DelegatedData(ex_policy))
{
}

inline void AphiGradDivAPenaltyCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    out_a_real_[index_i] -= a_divergence_penalty_ * grad_div_a_real_[index_i];
    out_a_imag_[index_i] -= a_divergence_penalty_ * grad_div_a_imag_[index_i];
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_A_DIVERGENCE_PENALTY_CK_HPP
