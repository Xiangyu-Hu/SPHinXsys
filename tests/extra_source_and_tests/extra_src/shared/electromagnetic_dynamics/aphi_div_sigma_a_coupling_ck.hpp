#ifndef APHI_DIV_SIGMA_A_COUPLING_CK_HPP
#define APHI_DIV_SIGMA_A_COUPLING_CK_HPP

#include "electromagnetic_dynamics/aphi_div_sigma_a_coupling_ck.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... Parameters>
inline AphiDivSigmaACouplingCK<Inner<Parameters...>>::AphiDivSigmaACouplingCK(
    Inner<Parameters...> &inner_relation, Real omega, const AphiVariableNames &variable_names)
    : BaseInteraction(inner_relation), omega_(omega),
      dv_a_real_(this->particles_->template getVariableByName<Vecd>(variable_names.solution.a_real)),
      dv_a_imag_(this->particles_->template getVariableByName<Vecd>(variable_names.solution.a_imag)),
      dv_lhs_phi_real_(this->particles_->template getVariableByName<Real>(variable_names.lhs.phi_real)),
      dv_lhs_phi_imag_(this->particles_->template getVariableByName<Real>(variable_names.lhs.phi_imag)),
      dv_sigma_(this->particles_->template getVariableByName<Real>(variable_names.material.sigma))
{
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiDivSigmaACouplingCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
      a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
      lhs_phi_real_(encloser.dv_lhs_phi_real_->DelegatedData(ex_policy)),
      lhs_phi_imag_(encloser.dv_lhs_phi_imag_->DelegatedData(ex_policy)),
      omega_(encloser.omega_)
{
}

template <typename... Parameters>
inline void AphiDivSigmaACouplingCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Vecd a_re_i = a_real_[index_i];
    const Vecd a_im_i = a_imag_[index_i];
    Real contribution_phi_real = 0.0;
    Real contribution_phi_imag = 0.0;

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        const Real sigma_ij = AphiHarmonicMean(sigma_i, sigma_[index_j]);
        contribution_phi_real += omega_ * sigma_ij * g_ij.dot(a_im_i - a_imag_[index_j]);
        contribution_phi_imag -= omega_ * sigma_ij * g_ij.dot(a_re_i - a_real_[index_j]);
    }

    lhs_phi_real_[index_i] += contribution_phi_real;
    lhs_phi_imag_[index_i] += contribution_phi_imag;
}

inline AphiDivSigmaAConstSigmaDiagnosticCouplingCK::AphiDivSigmaAConstSigmaDiagnosticCouplingCK(
    SPHBody &sph_body, Real omega, const AphiVariableNames &variable_names,
    const AphiDiagnosticNames &diagnostic_names)
    : LocalDynamics(sph_body), omega_(omega),
      dv_sigma_(particles_->template getVariableByName<Real>(variable_names.material.sigma)),
      dv_div_a_real_(particles_->template getVariableByName<Real>(diagnostic_names.div_a_real)),
      dv_div_a_imag_(particles_->template getVariableByName<Real>(diagnostic_names.div_a_imag)),
      dv_lhs_phi_real_(particles_->template getVariableByName<Real>(variable_names.lhs.phi_real)),
      dv_lhs_phi_imag_(particles_->template getVariableByName<Real>(variable_names.lhs.phi_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AphiDivSigmaAConstSigmaDiagnosticCouplingCK::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      div_a_real_(encloser.dv_div_a_real_->DelegatedData(ex_policy)),
      div_a_imag_(encloser.dv_div_a_imag_->DelegatedData(ex_policy)),
      lhs_phi_real_(encloser.dv_lhs_phi_real_->DelegatedData(ex_policy)),
      lhs_phi_imag_(encloser.dv_lhs_phi_imag_->DelegatedData(ex_policy)),
      omega_(encloser.omega_)
{
}

inline void AphiDivSigmaAConstSigmaDiagnosticCouplingCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    lhs_phi_real_[index_i] += omega_ * sigma_i * div_a_imag_[index_i];
    lhs_phi_imag_[index_i] -= omega_ * sigma_i * div_a_real_[index_i];
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_DIV_SIGMA_A_COUPLING_CK_HPP
