#ifndef APHI_JOULE_HEATING_CK_HPP
#define APHI_JOULE_HEATING_CK_HPP

#include "electromagnetic_dynamics/aphi_joule_heating_ck.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"

namespace SPH
{
namespace electromagnetics
{

inline RegisterAphiJouleHeatingFieldsCK::RegisterAphiJouleHeatingFieldsCK(SPHBody &sph_body,
                                                                          const AphiJouleHeatingFieldNames &field_names)
    : LocalDynamics(sph_body)
{
    auto &particles = particles_;
    particles->template registerStateVariable<Vecd>(field_names.grad_phi_real, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Vecd>(field_names.grad_phi_imag, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Vecd>(field_names.electric_field_a_real, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Vecd>(field_names.electric_field_a_imag, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Vecd>(field_names.current_density_real, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Vecd>(field_names.current_density_imag, ZeroData<Vecd>::value);
    particles->template registerStateVariable<Real>(field_names.joule_heat_source, Real(0));
    particles->addVariableToWrite<Real>(field_names.joule_heat_source);
}

template <typename... Parameters>
inline AphiComputeScalarPhiGradientCK<Inner<Parameters...>>::AphiComputeScalarPhiGradientCK(
    Inner<Parameters...> &inner_relation, const AphiBlockNames &phi_block, const AphiJouleHeatingFieldNames &field_names)
    : BaseInteraction(inner_relation)
{
    auto &particles = this->particles_;
    dv_phi_real_ = particles->template getVariableByName<Real>(phi_block.phi_real);
    dv_phi_imag_ = particles->template getVariableByName<Real>(phi_block.phi_imag);
    dv_grad_phi_real_ = particles->template getVariableByName<Vecd>(field_names.grad_phi_real);
    dv_grad_phi_imag_ = particles->template getVariableByName<Vecd>(field_names.grad_phi_imag);
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiComputeScalarPhiGradientCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
      grad_phi_real_(encloser.dv_grad_phi_real_->DelegatedData(ex_policy)),
      grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiComputeScalarPhiGradientCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real phi_re_i = phi_real_[index_i];
    const Real phi_im_i = phi_imag_[index_i];
    Vecd grad_phi_re(Real(0), Real(0), Real(0));
    Vecd grad_phi_im(Real(0), Real(0), Real(0));

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        grad_phi_re += g_ij * (phi_re_i - phi_real_[index_j]);
        grad_phi_im += g_ij * (phi_im_i - phi_imag_[index_j]);
    }

    grad_phi_real_[index_i] = grad_phi_re;
    grad_phi_imag_[index_i] = grad_phi_im;
}

template <typename... Parameters>
inline AphiComputeScalarPhiGradientCK<Contact<Parameters...>>::AphiComputeScalarPhiGradientCK(
    Contact<Parameters...> &contact_relation, const AphiBlockNames &phi_block, const AphiJouleHeatingFieldNames &field_names)
    : BaseInteraction(contact_relation)
{
    auto &particles = this->particles_;
    dv_phi_real_ = particles->template getVariableByName<Real>(phi_block.phi_real);
    dv_phi_imag_ = particles->template getVariableByName<Real>(phi_block.phi_imag);
    dv_grad_phi_real_ = particles->template getVariableByName<Vecd>(field_names.grad_phi_real);
    dv_grad_phi_imag_ = particles->template getVariableByName<Vecd>(field_names.grad_phi_imag);
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_phi_real_.push_back(contact_particles->template getVariableByName<Real>(phi_block.phi_real));
        dv_contact_phi_imag_.push_back(contact_particles->template getVariableByName<Real>(phi_block.phi_imag));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiComputeScalarPhiGradientCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
      contact_phi_real_(encloser.dv_contact_phi_real_[contact_index]->DelegatedData(ex_policy)),
      contact_phi_imag_(encloser.dv_contact_phi_imag_[contact_index]->DelegatedData(ex_policy)),
      grad_phi_real_(encloser.dv_grad_phi_real_->DelegatedData(ex_policy)),
      grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiComputeScalarPhiGradientCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real phi_re_i = phi_real_[index_i];
    const Real phi_im_i = phi_imag_[index_i];
    Vecd grad_phi_re = Vecd::Zero();
    Vecd grad_phi_im = Vecd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        grad_phi_re += g_ij * (phi_re_i - contact_phi_real_[index_j]);
        grad_phi_im += g_ij * (phi_im_i - contact_phi_imag_[index_j]);
    }

    grad_phi_real_[index_i] += grad_phi_re;
    grad_phi_imag_[index_i] += grad_phi_im;
}

inline AphiComputeFrequencyElectricFieldCK::AphiComputeFrequencyElectricFieldCK(
    SPHBody &sph_body, Real omega, const AphiBlockNames &solution_block, const AphiJouleHeatingFieldNames &field_names)
    : LocalDynamics(sph_body), omega_(omega)
{
    auto &particles = particles_;
    dv_a_real_ = particles->template getVariableByName<Vecd>(solution_block.a_real);
    dv_a_imag_ = particles->template getVariableByName<Vecd>(solution_block.a_imag);
    dv_grad_phi_real_ = particles->template getVariableByName<Vecd>(field_names.grad_phi_real);
    dv_grad_phi_imag_ = particles->template getVariableByName<Vecd>(field_names.grad_phi_imag);
    dv_electric_field_a_real_ = particles->template getVariableByName<Vecd>(field_names.electric_field_a_real);
    dv_electric_field_a_imag_ = particles->template getVariableByName<Vecd>(field_names.electric_field_a_imag);
}

template <class ExecutionPolicy, class EncloserType>
inline AphiComputeFrequencyElectricFieldCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                       EncloserType &encloser)
    : omega_(encloser.omega_), a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
      a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      grad_phi_real_(encloser.dv_grad_phi_real_->DelegatedData(ex_policy)),
      grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy)),
      electric_field_a_real_(encloser.dv_electric_field_a_real_->DelegatedData(ex_policy)),
      electric_field_a_imag_(encloser.dv_electric_field_a_imag_->DelegatedData(ex_policy))
{
}

inline void AphiComputeFrequencyElectricFieldCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    electric_field_a_real_[index_i] = omega_ * a_imag_[index_i] - grad_phi_real_[index_i];
    electric_field_a_imag_[index_i] = -omega_ * a_real_[index_i] - grad_phi_imag_[index_i];
}

inline AphiComputeJouleHeatSourceCK::AphiComputeJouleHeatSourceCK(SPHBody &sph_body,
                                                                  const AphiMaterialNames &material_names,
                                                                  const AphiJouleHeatingFieldNames &field_names)
    : LocalDynamics(sph_body)
{
    auto &particles = particles_;
    dv_sigma_ = particles->template getVariableByName<Real>(material_names.sigma);
    dv_electric_field_a_real_ = particles->template getVariableByName<Vecd>(field_names.electric_field_a_real);
    dv_electric_field_a_imag_ = particles->template getVariableByName<Vecd>(field_names.electric_field_a_imag);
    dv_current_density_real_ = particles->template getVariableByName<Vecd>(field_names.current_density_real);
    dv_current_density_imag_ = particles->template getVariableByName<Vecd>(field_names.current_density_imag);
    dv_joule_heat_source_ = particles->template getVariableByName<Real>(field_names.joule_heat_source);
}

template <class ExecutionPolicy, class EncloserType>
inline AphiComputeJouleHeatSourceCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                 EncloserType &encloser)
    : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      electric_field_a_real_(encloser.dv_electric_field_a_real_->DelegatedData(ex_policy)),
      electric_field_a_imag_(encloser.dv_electric_field_a_imag_->DelegatedData(ex_policy)),
      current_density_real_(encloser.dv_current_density_real_->DelegatedData(ex_policy)),
      current_density_imag_(encloser.dv_current_density_imag_->DelegatedData(ex_policy)),
      joule_heat_source_(encloser.dv_joule_heat_source_->DelegatedData(ex_policy))
{
}

inline void AphiComputeJouleHeatSourceCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd e_real = electric_field_a_real_[index_i];
    const Vecd e_imag = electric_field_a_imag_[index_i];
    const Real sigma_i = sigma_[index_i];
    const Vecd j_real = sigma_i * e_real;
    const Vecd j_imag = sigma_i * e_imag;
    current_density_real_[index_i] = j_real;
    current_density_imag_[index_i] = j_imag;
    joule_heat_source_[index_i] = Real(0.5) * (j_real.dot(e_real) + j_imag.dot(e_imag));
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_JOULE_HEATING_CK_HPP
