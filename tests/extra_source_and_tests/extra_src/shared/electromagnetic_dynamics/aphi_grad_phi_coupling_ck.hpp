#ifndef APHI_GRAD_PHI_COUPLING_CK_HPP
#define APHI_GRAD_PHI_COUPLING_CK_HPP

#include "electromagnetic_dynamics/aphi_grad_phi_coupling_ck.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... Parameters>
inline AphiGradPhiCouplingCK<Inner<Parameters...>>::AphiGradPhiCouplingCK(
    Inner<Parameters...> &inner_relation, const AphiVariableNames &variable_names)
    : BaseInteraction(inner_relation),
      dv_phi_real_(this->particles_->template getVariableByName<Real>(variable_names.solution.phi_real)),
      dv_phi_imag_(this->particles_->template getVariableByName<Real>(variable_names.solution.phi_imag)),
      dv_lhs_a_real_(this->particles_->template getVariableByName<Vecd>(variable_names.lhs.a_real)),
      dv_lhs_a_imag_(this->particles_->template getVariableByName<Vecd>(variable_names.lhs.a_imag)),
      dv_sigma_(this->particles_->template getVariableByName<Real>(variable_names.material.sigma))
{
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiGradPhiCouplingCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
      lhs_a_real_(encloser.dv_lhs_a_real_->DelegatedData(ex_policy)),
      lhs_a_imag_(encloser.dv_lhs_a_imag_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiGradPhiCouplingCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Real phi_re_i = phi_real_[index_i];
    const Real phi_im_i = phi_imag_[index_i];
    Vecd contribution_a_real = Vecd::Zero();
    Vecd contribution_a_imag = Vecd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        contribution_a_real += sigma_i * g_ij * (phi_re_i - phi_real_[index_j]);
        contribution_a_imag += sigma_i * g_ij * (phi_im_i - phi_imag_[index_j]);
    }

    lhs_a_real_[index_i] += contribution_a_real;
    lhs_a_imag_[index_i] += contribution_a_imag;
}

template <typename... Parameters>
inline AphiGradPhiCouplingCK<Contact<Parameters...>>::AphiGradPhiCouplingCK(
    Contact<Parameters...> &contact_relation, const AphiVariableNames &variable_names)
    : BaseInteraction(contact_relation),
      dv_phi_real_(this->particles_->template getVariableByName<Real>(variable_names.solution.phi_real)),
      dv_phi_imag_(this->particles_->template getVariableByName<Real>(variable_names.solution.phi_imag)),
      dv_lhs_a_real_(this->particles_->template getVariableByName<Vecd>(variable_names.lhs.a_real)),
      dv_lhs_a_imag_(this->particles_->template getVariableByName<Vecd>(variable_names.lhs.a_imag)),
      dv_sigma_(this->particles_->template getVariableByName<Real>(variable_names.material.sigma))
{
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_phi_real_.push_back(contact_particles->template getVariableByName<Real>(variable_names.solution.phi_real));
        dv_contact_phi_imag_.push_back(contact_particles->template getVariableByName<Real>(variable_names.solution.phi_imag));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiGradPhiCouplingCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
      phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
      contact_phi_real_(encloser.dv_contact_phi_real_[contact_index]->DelegatedData(ex_policy)),
      contact_phi_imag_(encloser.dv_contact_phi_imag_[contact_index]->DelegatedData(ex_policy)),
      lhs_a_real_(encloser.dv_lhs_a_real_->DelegatedData(ex_policy)),
      lhs_a_imag_(encloser.dv_lhs_a_imag_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiGradPhiCouplingCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Real phi_re_i = phi_real_[index_i];
    const Real phi_im_i = phi_imag_[index_i];
    Vecd contribution_a_real = Vecd::Zero();
    Vecd contribution_a_imag = Vecd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        contribution_a_real += sigma_i * g_ij * (phi_re_i - contact_phi_real_[index_j]);
        contribution_a_imag += sigma_i * g_ij * (phi_im_i - contact_phi_imag_[index_j]);
    }

    lhs_a_real_[index_i] += contribution_a_real;
    lhs_a_imag_[index_i] += contribution_a_imag;
}

inline AphiGradPhiDiagnosticCouplingCK::AphiGradPhiDiagnosticCouplingCK(SPHBody &sph_body,
                                                                        const AphiVariableNames &variable_names)
    : LocalDynamics(sph_body),
      dv_sigma_(particles_->template getVariableByName<Real>(variable_names.material.sigma)),
      dv_grad_phi_real_(particles_->template getVariableByName<Vecd>(variable_names.solution.phi_real + "Gradient")),
      dv_grad_phi_imag_(particles_->template getVariableByName<Vecd>(variable_names.solution.phi_imag + "Gradient")),
      dv_lhs_a_real_(particles_->template getVariableByName<Vecd>(variable_names.lhs.a_real)),
      dv_lhs_a_imag_(particles_->template getVariableByName<Vecd>(variable_names.lhs.a_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiGradPhiDiagnosticCouplingCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                                                                   EncloserType &encloser)
    : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      grad_phi_real_(encloser.dv_grad_phi_real_->DelegatedData(ex_policy)),
      grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy)),
      lhs_a_real_(encloser.dv_lhs_a_real_->DelegatedData(ex_policy)),
      lhs_a_imag_(encloser.dv_lhs_a_imag_->DelegatedData(ex_policy))
{
}

inline void AphiGradPhiDiagnosticCouplingCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    lhs_a_real_[index_i] += sigma_i * grad_phi_real_[index_i];
    lhs_a_imag_[index_i] += sigma_i * grad_phi_imag_[index_i];
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GRAD_PHI_COUPLING_CK_HPP
