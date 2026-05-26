#ifndef APHI_REACTION_CK_HPP
#define APHI_REACTION_CK_HPP

#include "electromagnetic_dynamics/aphi_reaction_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline AphiReactionCK::AphiReactionCK(SPHBody &sph_body, Real omega, const AphiVariableNames &variable_names)
    : LocalDynamics(sph_body), omega_(omega),
      dv_input_a_real_(particles_->template getVariableByName<Vecd>(variable_names.solution.a_real)),
      dv_input_a_imag_(particles_->template getVariableByName<Vecd>(variable_names.solution.a_imag)),
      dv_output_a_real_(particles_->template getVariableByName<Vecd>(variable_names.lhs.a_real)),
      dv_output_a_imag_(particles_->template getVariableByName<Vecd>(variable_names.lhs.a_imag)),
      dv_sigma_(particles_->template getVariableByName<Real>(variable_names.material.sigma))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiReactionCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : input_a_real_(encloser.dv_input_a_real_->DelegatedData(ex_policy)),
      input_a_imag_(encloser.dv_input_a_imag_->DelegatedData(ex_policy)),
      output_a_real_(encloser.dv_output_a_real_->DelegatedData(ex_policy)),
      output_a_imag_(encloser.dv_output_a_imag_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      omega_(encloser.omega_)
{
}

inline void AphiReactionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    output_a_real_[index_i] += -omega_ * sigma_i * input_a_imag_[index_i];
    output_a_imag_[index_i] += omega_ * sigma_i * input_a_real_[index_i];
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_REACTION_CK_HPP
