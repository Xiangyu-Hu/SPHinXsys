#ifndef ELECTROMAGNETIC_OPHELIE_POSTPROCESS_HPP
#define ELECTROMAGNETIC_OPHELIE_POSTPROCESS_HPP

#include "electromagnetic_ophelie_postprocess.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline ComputeOphelieEJQFromASrcNoPhiCK::ComputeOphelieEJQFromASrcNoPhiCK(SPHBody &sph_body,
                                                                            const OphelieGlassFieldNames &names,
                                                                            const OphelieParameters &params)
    : LocalDynamics(sph_body), omega_(params.omega()),
      dv_sigma_(particles_->template getVariableByName<Real>(names.sigma)),
      dv_a_src_real_(particles_->template getVariableByName<Vecd>(names.a_src_real)),
      dv_a_src_imag_(particles_->template getVariableByName<Vecd>(names.a_src_imag)),
      dv_e_real_(particles_->template getVariableByName<Vecd>(names.e_real)),
      dv_e_imag_(particles_->template getVariableByName<Vecd>(names.e_imag)),
      dv_j_real_(particles_->template getVariableByName<Vecd>(names.j_real)),
      dv_j_imag_(particles_->template getVariableByName<Vecd>(names.j_imag)),
      dv_joule_heat_(particles_->template getVariableByName<Real>(names.joule_heat))
{
}

template <class ExecutionPolicy, class EncloserType>
inline ComputeOphelieEJQFromASrcNoPhiCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                    EncloserType &encloser)
    : omega_(encloser.omega_), sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      a_src_real_(encloser.dv_a_src_real_->DelegatedData(ex_policy)),
      a_src_imag_(encloser.dv_a_src_imag_->DelegatedData(ex_policy)),
      e_real_(encloser.dv_e_real_->DelegatedData(ex_policy)),
      e_imag_(encloser.dv_e_imag_->DelegatedData(ex_policy)),
      j_real_(encloser.dv_j_real_->DelegatedData(ex_policy)),
      j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)),
      joule_heat_(encloser.dv_joule_heat_->DelegatedData(ex_policy))
{
}

inline void ComputeOphelieEJQFromASrcNoPhiCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Vecd e_real = omega_ * a_src_imag_[index_i];
    const Vecd e_imag = -omega_ * a_src_real_[index_i];
    const Vecd j_real = sigma_i * e_real;
    const Vecd j_imag = sigma_i * e_imag;

    e_real_[index_i] = e_real;
    e_imag_[index_i] = e_imag;
    j_real_[index_i] = j_real;
    j_imag_[index_i] = j_imag;
    joule_heat_[index_i] = Real(0.5) * (j_real.dot(e_real) + j_imag.dot(e_imag));
}

inline ScaleOphelieElectromagneticFieldsCK::ScaleOphelieElectromagneticFieldsCK(SPHBody &sph_body,
                                                                               const OphelieGlassFieldNames &names,
                                                                               Real field_scale, Real power_scale)
    : LocalDynamics(sph_body), field_scale_(field_scale), power_scale_(power_scale),
      dv_a_coil_real_(particles_->template getVariableByName<Vecd>(names.a_coil_real)),
      dv_a_coil_imag_(particles_->template getVariableByName<Vecd>(names.a_coil_imag)),
      dv_b_coil_real_(particles_->template getVariableByName<Vecd>(names.b_coil_real)),
      dv_b_coil_imag_(particles_->template getVariableByName<Vecd>(names.b_coil_imag)),
      dv_a_ind_real_(particles_->template getVariableByName<Vecd>(names.a_ind_real)),
      dv_a_ind_imag_(particles_->template getVariableByName<Vecd>(names.a_ind_imag)),
      dv_b_ind_real_(particles_->template getVariableByName<Vecd>(names.b_ind_real)),
      dv_b_ind_imag_(particles_->template getVariableByName<Vecd>(names.b_ind_imag)),
      dv_a_src_real_(particles_->template getVariableByName<Vecd>(names.a_src_real)),
      dv_a_src_imag_(particles_->template getVariableByName<Vecd>(names.a_src_imag)),
      dv_b_src_real_(particles_->template getVariableByName<Vecd>(names.b_src_real)),
      dv_b_src_imag_(particles_->template getVariableByName<Vecd>(names.b_src_imag)),
      dv_e_real_(particles_->template getVariableByName<Vecd>(names.e_real)),
      dv_e_imag_(particles_->template getVariableByName<Vecd>(names.e_imag)),
      dv_j_real_(particles_->template getVariableByName<Vecd>(names.j_real)),
      dv_j_imag_(particles_->template getVariableByName<Vecd>(names.j_imag)),
      dv_phi_imag_(particles_->template getVariableByName<Real>(names.phi_imag)),
      dv_grad_phi_imag_(particles_->template getVariableByName<Vecd>(names.grad_phi_imag)),
      dv_joule_heat_(particles_->template getVariableByName<Real>(names.joule_heat))
{
}

template <class ExecutionPolicy, class EncloserType>
inline ScaleOphelieElectromagneticFieldsCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                       EncloserType &encloser)
    : field_scale_(encloser.field_scale_), power_scale_(encloser.power_scale_),
      a_coil_real_(encloser.dv_a_coil_real_->DelegatedData(ex_policy)),
      a_coil_imag_(encloser.dv_a_coil_imag_->DelegatedData(ex_policy)),
      b_coil_real_(encloser.dv_b_coil_real_->DelegatedData(ex_policy)),
      b_coil_imag_(encloser.dv_b_coil_imag_->DelegatedData(ex_policy)),
      a_ind_real_(encloser.dv_a_ind_real_->DelegatedData(ex_policy)),
      a_ind_imag_(encloser.dv_a_ind_imag_->DelegatedData(ex_policy)),
      b_ind_real_(encloser.dv_b_ind_real_->DelegatedData(ex_policy)),
      b_ind_imag_(encloser.dv_b_ind_imag_->DelegatedData(ex_policy)),
      a_src_real_(encloser.dv_a_src_real_->DelegatedData(ex_policy)),
      a_src_imag_(encloser.dv_a_src_imag_->DelegatedData(ex_policy)),
      b_src_real_(encloser.dv_b_src_real_->DelegatedData(ex_policy)),
      b_src_imag_(encloser.dv_b_src_imag_->DelegatedData(ex_policy)),
      e_real_(encloser.dv_e_real_->DelegatedData(ex_policy)),
      e_imag_(encloser.dv_e_imag_->DelegatedData(ex_policy)),
      j_real_(encloser.dv_j_real_->DelegatedData(ex_policy)),
      j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
      grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy)),
      joule_heat_(encloser.dv_joule_heat_->DelegatedData(ex_policy))
{
}

inline void ScaleOphelieElectromagneticFieldsCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const auto scale_vec = [this](Vecd &value) { value *= field_scale_; };
    scale_vec(a_coil_real_[index_i]);
    scale_vec(a_coil_imag_[index_i]);
    scale_vec(b_coil_real_[index_i]);
    scale_vec(b_coil_imag_[index_i]);
    scale_vec(a_ind_real_[index_i]);
    scale_vec(a_ind_imag_[index_i]);
    scale_vec(b_ind_real_[index_i]);
    scale_vec(b_ind_imag_[index_i]);
    scale_vec(a_src_real_[index_i]);
    scale_vec(a_src_imag_[index_i]);
    scale_vec(b_src_real_[index_i]);
    scale_vec(b_src_imag_[index_i]);
    scale_vec(e_real_[index_i]);
    scale_vec(e_imag_[index_i]);
    scale_vec(j_real_[index_i]);
    scale_vec(j_imag_[index_i]);
    grad_phi_imag_[index_i] *= field_scale_;
    phi_imag_[index_i] *= field_scale_;
    joule_heat_[index_i] *= power_scale_;
}

inline ScaleOphelieJouleHeatCK::ScaleOphelieJouleHeatCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, Real scale)
    : LocalDynamics(sph_body), scale_(scale), dv_joule_heat_(particles_->template getVariableByName<Real>(names.joule_heat))
{
}

template <class ExecutionPolicy, class EncloserType>
inline ScaleOphelieJouleHeatCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : scale_(encloser.scale_), joule_heat_(encloser.dv_joule_heat_->DelegatedData(ex_policy))
{
}

inline void ScaleOphelieJouleHeatCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    joule_heat_[index_i] *= scale_;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_POSTPROCESS_HPP
