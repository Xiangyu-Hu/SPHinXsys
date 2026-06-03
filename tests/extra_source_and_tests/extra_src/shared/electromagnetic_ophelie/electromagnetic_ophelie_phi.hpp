#ifndef ELECTROMAGNETIC_OPHELIE_PHI_HPP
#define ELECTROMAGNETIC_OPHELIE_PHI_HPP

#include "electromagnetic_ophelie_phi.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline ComputeOphelieEJQWithPhiCK::ComputeOphelieEJQWithPhiCK(SPHBody &sph_body, const OphelieGlassFieldNames &names,
                                                              const OphelieParameters &params)
    : LocalDynamics(sph_body), omega_(params.omega()),
      dv_sigma_(particles_->template getVariableByName<Real>(names.sigma)),
      dv_a_src_real_(particles_->template getVariableByName<Vecd>(names.a_src_real)),
      dv_grad_phi_imag_(particles_->template getVariableByName<Vecd>(names.grad_phi_imag)),
      dv_e_real_(particles_->template getVariableByName<Vecd>(names.e_real)),
      dv_e_imag_(particles_->template getVariableByName<Vecd>(names.e_imag)),
      dv_j_real_(particles_->template getVariableByName<Vecd>(names.j_real)),
      dv_j_imag_(particles_->template getVariableByName<Vecd>(names.j_imag)),
      dv_joule_heat_(particles_->template getVariableByName<Real>(names.joule_heat))
{
}

template <class ExecutionPolicy, class EncloserType>
inline ComputeOphelieEJQWithPhiCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : omega_(encloser.omega_), sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      a_src_real_(encloser.dv_a_src_real_->DelegatedData(ex_policy)),
      grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy)),
      e_real_(encloser.dv_e_real_->DelegatedData(ex_policy)), e_imag_(encloser.dv_e_imag_->DelegatedData(ex_policy)),
      j_real_(encloser.dv_j_real_->DelegatedData(ex_policy)), j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)),
      joule_heat_(encloser.dv_joule_heat_->DelegatedData(ex_policy))
{
}

inline void ComputeOphelieEJQWithPhiCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Vecd e_real = Vecd::Zero();
    const Vecd e_imag = -grad_phi_imag_[index_i] - omega_ * a_src_real_[index_i];
    const Vecd j_real = sigma_i * e_real;
    const Vecd j_imag = sigma_i * e_imag;
    e_real_[index_i] = e_real;
    e_imag_[index_i] = e_imag;
    j_real_[index_i] = j_real;
    j_imag_[index_i] = j_imag;
    joule_heat_[index_i] = Real(0.5) * (j_real.dot(e_real) + j_imag.dot(e_imag));
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_HPP
