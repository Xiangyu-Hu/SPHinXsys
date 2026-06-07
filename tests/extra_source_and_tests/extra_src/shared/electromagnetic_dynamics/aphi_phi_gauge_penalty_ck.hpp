#ifndef APHI_PHI_GAUGE_PENALTY_CK_HPP
#define APHI_PHI_GAUGE_PENALTY_CK_HPP

#include "electromagnetic_dynamics/aphi_phi_gauge_penalty_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiPhiGaugePenaltyCK::AphiPhiGaugePenaltyCK(SPHBody &sph_body, const AphiBlockNames &input_block,
                                                  const AphiBlockNames &output_block, Real phi_gauge_penalty)
    : LocalDynamics(sph_body), phi_gauge_penalty_(phi_gauge_penalty),
      dv_in_phi_real_(particles_->template getVariableByName<Real>(input_block.phi_real)),
      dv_in_phi_imag_(particles_->template getVariableByName<Real>(input_block.phi_imag)),
      dv_out_phi_real_(particles_->template getVariableByName<Real>(output_block.phi_real)),
      dv_out_phi_imag_(particles_->template getVariableByName<Real>(output_block.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiPhiGaugePenaltyCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : phi_gauge_penalty_(encloser.phi_gauge_penalty_),
      in_phi_real_(encloser.dv_in_phi_real_->DelegatedData(ex_policy)),
      in_phi_imag_(encloser.dv_in_phi_imag_->DelegatedData(ex_policy)),
      out_phi_real_(encloser.dv_out_phi_real_->DelegatedData(ex_policy)),
      out_phi_imag_(encloser.dv_out_phi_imag_->DelegatedData(ex_policy))
{
}

inline void AphiPhiGaugePenaltyCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    out_phi_real_[index_i] += phi_gauge_penalty_ * in_phi_real_[index_i];
    out_phi_imag_[index_i] += phi_gauge_penalty_ * in_phi_imag_[index_i];
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PHI_GAUGE_PENALTY_CK_HPP
