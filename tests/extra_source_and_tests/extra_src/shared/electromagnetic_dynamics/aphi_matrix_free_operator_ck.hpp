#ifndef APHI_MATRIX_FREE_OPERATOR_CK_HPP
#define APHI_MATRIX_FREE_OPERATOR_CK_HPP

#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"

#include <stdexcept>

namespace SPH
{
namespace electromagnetics
{

template <typename... Parameters>
inline AphiApplyCK<Inner<Parameters...>>::AphiApplyCK(
    Inner<Parameters...> &inner_relation, const AphiBlockNames &input_block, const AphiBlockNames &output_block,
    const AphiMaterialNames &material_names, Real omega, const AphiLhsAssemblyOptions &options,
    Real pair_weight_regularization)
    : BaseInteraction(inner_relation),
      options_(options),
      pair_weight_regularization_(pair_weight_regularization),
      reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      dv_in_a_real_(this->particles_->template getVariableByName<Vecd>(input_block.a_real)),
      dv_in_a_imag_(this->particles_->template getVariableByName<Vecd>(input_block.a_imag)),
      dv_in_phi_real_(this->particles_->template getVariableByName<Real>(input_block.phi_real)),
      dv_in_phi_imag_(this->particles_->template getVariableByName<Real>(input_block.phi_imag)),
      dv_out_a_real_(this->particles_->template getVariableByName<Vecd>(output_block.a_real)),
      dv_out_a_imag_(this->particles_->template getVariableByName<Vecd>(output_block.a_imag)),
      dv_out_phi_real_(this->particles_->template getVariableByName<Real>(output_block.phi_real)),
      dv_out_phi_imag_(this->particles_->template getVariableByName<Real>(output_block.phi_imag)),
      dv_sigma_(this->particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(this->particles_->template getVariableByName<Real>(material_names.nu))
{
    options_.omega = omega;
    if (options_.grad_phi_mode != AphiGradPhiCouplingMode::PairwiseUncorrected ||
        options_.div_mode != AphiDivCouplingMode::PairwiseUncorrected)
    {
        throw std::runtime_error(
            "AphiApplyCK supports only PairwiseUncorrected grad/div modes. "
            "Use AphiAssembleLhsDebugDynamicsBundle for B-corrected diagnostic paths.");
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiApplyCK<Inner<Parameters...>>::InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy,
                                                                         EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
      nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      in_a_real_(encloser.dv_in_a_real_->DelegatedData(ex_policy)),
      in_a_imag_(encloser.dv_in_a_imag_->DelegatedData(ex_policy)),
      in_phi_real_(encloser.dv_in_phi_real_->DelegatedData(ex_policy)),
      in_phi_imag_(encloser.dv_in_phi_imag_->DelegatedData(ex_policy)),
      out_a_real_(encloser.dv_out_a_real_->DelegatedData(ex_policy)),
      out_a_imag_(encloser.dv_out_a_imag_->DelegatedData(ex_policy)),
      out_phi_real_(encloser.dv_out_phi_real_->DelegatedData(ex_policy)),
      out_phi_imag_(encloser.dv_out_phi_imag_->DelegatedData(ex_policy)),
      terms_(encloser.options_.terms),
      omega_(encloser.options_.omega),
      use_phi_gauge_penalty_(encloser.options_.use_phi_gauge_penalty),
      phi_gauge_penalty_(encloser.options_.phi_gauge_penalty),
      pair_weight_regularization_(encloser.pair_weight_regularization_),
      reference_smoothing_length_(encloser.reference_smoothing_length_)
{
}

template <typename... Parameters>
inline void AphiApplyCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Real nu_i = nu_[index_i];
    const Vecd a_re_i = in_a_real_[index_i];
    const Vecd a_im_i = in_a_imag_[index_i];
    const Real phi_re_i = in_phi_real_[index_i];
    const Real phi_im_i = in_phi_imag_[index_i];

    Vecd lhs_a_real = Vecd::Zero();
    Vecd lhs_a_imag = Vecd::Zero();
    Real lhs_phi_real = Real(0);
    Real lhs_phi_imag = Real(0);

    if (terms_.laplace_phi)
    {
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
            const Real distance = r_ij_vec.norm();
            const Real distance_sq = r_ij_vec.squaredNorm();
            const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            const Real sigma_ij = AphiHarmonicMean(sigma_i, sigma_[index_j]);
            const Real laplace_weight =
                sigma_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                             reference_smoothing_length_);
            lhs_phi_real += laplace_weight * (phi_re_i - in_phi_real_[index_j]);
            lhs_phi_imag += laplace_weight * (phi_im_i - in_phi_imag_[index_j]);
        }
    }

    if (terms_.laplace_a)
    {
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
            const Real distance = r_ij_vec.norm();
            const Real distance_sq = r_ij_vec.squaredNorm();
            const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            const Real nu_ij = AphiHarmonicMean(nu_i, nu_[index_j]);
            const Real laplace_weight =
                nu_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                          reference_smoothing_length_);
            lhs_a_real += laplace_weight * (a_re_i - in_a_real_[index_j]);
            lhs_a_imag += laplace_weight * (a_im_i - in_a_imag_[index_j]);
        }
    }

    if (terms_.reaction)
    {
        lhs_a_real += -omega_ * sigma_i * a_im_i;
        lhs_a_imag += omega_ * sigma_i * a_re_i;
    }

    if (terms_.grad_phi_coupling)
    {
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
            lhs_a_real += sigma_i * g_ij * (phi_re_i - in_phi_real_[index_j]);
            lhs_a_imag += sigma_i * g_ij * (phi_im_i - in_phi_imag_[index_j]);
        }
    }

    if (terms_.div_sigma_a_coupling)
    {
        for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        {
            const UnsignedInt index_j = this->neighbor_index_[n];
            const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
            const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
            const Real sigma_ij = AphiHarmonicMean(sigma_i, sigma_[index_j]);
            lhs_phi_real += omega_ * sigma_ij * g_ij.dot(a_im_i - in_a_imag_[index_j]);
            lhs_phi_imag -= omega_ * sigma_ij * g_ij.dot(a_re_i - in_a_real_[index_j]);
        }
    }

    if (use_phi_gauge_penalty_)
    {
        lhs_phi_real += phi_gauge_penalty_ * phi_re_i;
        lhs_phi_imag += phi_gauge_penalty_ * phi_im_i;
    }

    out_a_real_[index_i] = lhs_a_real;
    out_a_imag_[index_i] = lhs_a_imag;
    out_phi_real_[index_i] = lhs_phi_real;
    out_phi_imag_[index_i] = lhs_phi_imag;
}

template <typename... Parameters>
inline AphiApplyCK<Contact<Parameters...>>::AphiApplyCK(
    Contact<Parameters...> &contact_relation, const AphiBlockNames &input_block, const AphiBlockNames &output_block,
    const AphiMaterialNames &material_names, Real omega, const AphiLhsAssemblyOptions &options,
    Real pair_weight_regularization)
    : BaseInteraction(contact_relation),
      options_(options),
      pair_weight_regularization_(pair_weight_regularization),
      reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      dv_in_a_real_(this->particles_->template getVariableByName<Vecd>(input_block.a_real)),
      dv_in_a_imag_(this->particles_->template getVariableByName<Vecd>(input_block.a_imag)),
      dv_in_phi_real_(this->particles_->template getVariableByName<Real>(input_block.phi_real)),
      dv_in_phi_imag_(this->particles_->template getVariableByName<Real>(input_block.phi_imag)),
      dv_out_a_real_(this->particles_->template getVariableByName<Vecd>(output_block.a_real)),
      dv_out_a_imag_(this->particles_->template getVariableByName<Vecd>(output_block.a_imag)),
      dv_out_phi_real_(this->particles_->template getVariableByName<Real>(output_block.phi_real)),
      dv_out_phi_imag_(this->particles_->template getVariableByName<Real>(output_block.phi_imag)),
      dv_sigma_(this->particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(this->particles_->template getVariableByName<Real>(material_names.nu))
{
    options_.omega = omega;
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_sigma_.push_back(contact_particles->template getVariableByName<Real>(material_names.sigma));
        dv_contact_nu_.push_back(contact_particles->template getVariableByName<Real>(material_names.nu));
        dv_contact_in_phi_real_.push_back(contact_particles->template getVariableByName<Real>(input_block.phi_real));
        dv_contact_in_phi_imag_.push_back(contact_particles->template getVariableByName<Real>(input_block.phi_imag));
        dv_contact_in_a_real_.push_back(contact_particles->template getVariableByName<Vecd>(input_block.a_real));
        dv_contact_in_a_imag_.push_back(contact_particles->template getVariableByName<Vecd>(input_block.a_imag));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiApplyCK<Contact<Parameters...>>::InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy,
                                                                           EncloserType &encloser,
                                                                           size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      contact_sigma_(encloser.dv_contact_sigma_[contact_index]->DelegatedData(ex_policy)),
      contact_nu_(encloser.dv_contact_nu_[contact_index]->DelegatedData(ex_policy)),
      in_phi_real_(encloser.dv_in_phi_real_->DelegatedData(ex_policy)),
      in_phi_imag_(encloser.dv_in_phi_imag_->DelegatedData(ex_policy)),
      in_a_real_(encloser.dv_in_a_real_->DelegatedData(ex_policy)),
      in_a_imag_(encloser.dv_in_a_imag_->DelegatedData(ex_policy)),
      contact_in_phi_real_(encloser.dv_contact_in_phi_real_[contact_index]->DelegatedData(ex_policy)),
      contact_in_phi_imag_(encloser.dv_contact_in_phi_imag_[contact_index]->DelegatedData(ex_policy)),
      contact_in_a_real_(encloser.dv_contact_in_a_real_[contact_index]->DelegatedData(ex_policy)),
      contact_in_a_imag_(encloser.dv_contact_in_a_imag_[contact_index]->DelegatedData(ex_policy)),
      out_a_real_(encloser.dv_out_a_real_->DelegatedData(ex_policy)),
      out_a_imag_(encloser.dv_out_a_imag_->DelegatedData(ex_policy)),
      out_phi_real_(encloser.dv_out_phi_real_->DelegatedData(ex_policy)),
      out_phi_imag_(encloser.dv_out_phi_imag_->DelegatedData(ex_policy)), terms_(encloser.options_.terms),
      omega_(encloser.options_.omega), pair_weight_regularization_(encloser.pair_weight_regularization_),
      reference_smoothing_length_(encloser.reference_smoothing_length_)
{
}

template <typename... Parameters>
inline void AphiApplyCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Real nu_i = nu_[index_i];
    const Vecd a_re_i = in_a_real_[index_i];
    const Vecd a_im_i = in_a_imag_[index_i];
    const Real phi_re_i = in_phi_real_[index_i];
    const Real phi_im_i = in_phi_imag_[index_i];

    Vecd lhs_a_real = Vecd::Zero();
    Vecd lhs_a_imag = Vecd::Zero();
    Real lhs_phi_real = Real(0);
    Real lhs_phi_imag = Real(0);

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
        const Real distance = r_ij_vec.norm();
        const Real distance_sq = r_ij_vec.squaredNorm();
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];

        if (terms_.laplace_phi)
        {
            const Real sigma_ij = AphiHarmonicMean(sigma_i, contact_sigma_[index_j]);
            const Real laplace_weight =
                sigma_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                             reference_smoothing_length_);
            lhs_phi_real += laplace_weight * (phi_re_i - contact_in_phi_real_[index_j]);
            lhs_phi_imag += laplace_weight * (phi_im_i - contact_in_phi_imag_[index_j]);
        }

        if (terms_.laplace_a)
        {
            const Real nu_ij = AphiHarmonicMean(nu_i, contact_nu_[index_j]);
            const Real laplace_weight =
                nu_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                           reference_smoothing_length_);
            lhs_a_real += laplace_weight * (a_re_i - contact_in_a_real_[index_j]);
            lhs_a_imag += laplace_weight * (a_im_i - contact_in_a_imag_[index_j]);
        }

        if (terms_.grad_phi_coupling || terms_.div_sigma_a_coupling)
        {
            const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
            if (terms_.grad_phi_coupling)
            {
                lhs_a_real += sigma_i * g_ij * (phi_re_i - contact_in_phi_real_[index_j]);
                lhs_a_imag += sigma_i * g_ij * (phi_im_i - contact_in_phi_imag_[index_j]);
            }
            if (terms_.div_sigma_a_coupling)
            {
                const Real sigma_ij = AphiHarmonicMean(sigma_i, contact_sigma_[index_j]);
                lhs_phi_real += omega_ * sigma_ij * g_ij.dot(a_im_i - contact_in_a_imag_[index_j]);
                lhs_phi_imag -= omega_ * sigma_ij * g_ij.dot(a_re_i - contact_in_a_real_[index_j]);
            }
        }
    }

    out_a_real_[index_i] += lhs_a_real;
    out_a_imag_[index_i] += lhs_a_imag;
    out_phi_real_[index_i] += lhs_phi_real;
    out_phi_imag_[index_i] += lhs_phi_imag;
}

template <class ExecutionPolicy>
inline AphiApplyDynamicsBundle<ExecutionPolicy>::AphiApplyDynamicsBundle(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiBlockNames &input_block,
    const AphiBlockNames &output_block, const AphiMaterialNames &material_names, Real omega,
    const AphiLhsAssemblyOptions &options, Real pair_weight_regularization)
    : options_(options),
      input_block_(input_block),
      output_block_(output_block),
      sph_body_(sph_body),
      inner_relation_(inner_relation),
      zero_output_(sph_body, output_block),
      apply_(DynamicsArgs(inner_relation, input_block, output_block, material_names, omega, options,
                          pair_weight_regularization)),
      a_divergence_penalty_pairwise_(sph_body, inner_relation, input_block, output_block, options.a_divergence_penalty),
      a_divergence_penalty_b_corrected_(sph_body, inner_relation, input_block, output_block,
                                          options.a_divergence_penalty)
{
}

template <class ExecutionPolicy>
inline void AphiApplyDynamicsBundle<ExecutionPolicy>::exec()
{
    zero_output_.exec();
    apply_.exec();
    if (options_.use_a_divergence_penalty)
    {
        if (options_.a_divergence_penalty_mode == AphiADivergencePenaltyMode::PairwiseUncorrected)
        {
            a_divergence_penalty_pairwise_.exec();
        }
        else
        {
            a_divergence_penalty_b_corrected_.exec();
        }
    }
}

template <class ExecutionPolicy>
inline AphiApplyContactDynamicsBundle<ExecutionPolicy>::AphiApplyContactDynamicsBundle(
    SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation, const AphiBlockNames &input_block,
    const AphiBlockNames &output_block, const AphiMaterialNames &material_names, Real omega,
    const AphiLhsAssemblyOptions &options, Real pair_weight_regularization)
    : options_(options),
      zero_output_(sph_body, output_block),
      apply_inner_(DynamicsArgs(inner_relation, input_block, output_block, material_names, omega, options,
                              pair_weight_regularization)),
      apply_contact_(DynamicsArgs(contact_relation, input_block, output_block, material_names, omega, options,
                                pair_weight_regularization)),
      a_divergence_penalty_inner_contact_(sph_body, inner_relation, contact_relation, input_block, output_block,
                                          options.a_divergence_penalty),
      a_divergence_penalty_inner_only_(sph_body, inner_relation, input_block, output_block, options.a_divergence_penalty)
{
}

template <class ExecutionPolicy>
inline void AphiApplyContactDynamicsBundle<ExecutionPolicy>::exec()
{
    zero_output_.exec();
    apply_inner_.exec();
    apply_contact_.exec();
    if (options_.use_a_divergence_penalty &&
        options_.a_divergence_penalty_mode == AphiADivergencePenaltyMode::PairwiseUncorrected)
    {
        if (options_.contact_a_divergence_penalty_stencil == AphiContactADivergencePenaltyStencilMode::InnerOnly)
        {
            a_divergence_penalty_inner_only_.exec();
        }
        else
        {
            a_divergence_penalty_inner_contact_.exec();
        }
    }
}

template <typename... Parameters>
inline AphiApplyContactBlockDiagonalCK<Contact<Parameters...>>::AphiApplyContactBlockDiagonalCK(
    Contact<Parameters...> &contact_relation, const AphiBlockNames &input_block, const AphiBlockNames &output_block,
    const AphiMaterialNames &material_names, Real omega, const AphiLhsAssemblyOptions &options,
    Real pair_weight_regularization)
    : BaseInteraction(contact_relation),
      options_(options),
      pair_weight_regularization_(pair_weight_regularization),
      reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      dv_in_a_real_(this->particles_->template getVariableByName<Vecd>(input_block.a_real)),
      dv_in_a_imag_(this->particles_->template getVariableByName<Vecd>(input_block.a_imag)),
      dv_in_phi_real_(this->particles_->template getVariableByName<Real>(input_block.phi_real)),
      dv_in_phi_imag_(this->particles_->template getVariableByName<Real>(input_block.phi_imag)),
      dv_out_a_real_(this->particles_->template getVariableByName<Vecd>(output_block.a_real)),
      dv_out_a_imag_(this->particles_->template getVariableByName<Vecd>(output_block.a_imag)),
      dv_out_phi_real_(this->particles_->template getVariableByName<Real>(output_block.phi_real)),
      dv_out_phi_imag_(this->particles_->template getVariableByName<Real>(output_block.phi_imag)),
      dv_sigma_(this->particles_->template getVariableByName<Real>(material_names.sigma)),
      dv_nu_(this->particles_->template getVariableByName<Real>(material_names.nu))
{
    options_.omega = omega;
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_sigma_.push_back(contact_particles->template getVariableByName<Real>(material_names.sigma));
        dv_contact_nu_.push_back(contact_particles->template getVariableByName<Real>(material_names.nu));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiApplyContactBlockDiagonalCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), nu_(encloser.dv_nu_->DelegatedData(ex_policy)),
      contact_sigma_(encloser.dv_contact_sigma_[contact_index]->DelegatedData(ex_policy)),
      contact_nu_(encloser.dv_contact_nu_[contact_index]->DelegatedData(ex_policy)),
      in_phi_real_(encloser.dv_in_phi_real_->DelegatedData(ex_policy)),
      in_phi_imag_(encloser.dv_in_phi_imag_->DelegatedData(ex_policy)),
      in_a_real_(encloser.dv_in_a_real_->DelegatedData(ex_policy)),
      in_a_imag_(encloser.dv_in_a_imag_->DelegatedData(ex_policy)),
      out_a_real_(encloser.dv_out_a_real_->DelegatedData(ex_policy)),
      out_a_imag_(encloser.dv_out_a_imag_->DelegatedData(ex_policy)),
      out_phi_real_(encloser.dv_out_phi_real_->DelegatedData(ex_policy)),
      out_phi_imag_(encloser.dv_out_phi_imag_->DelegatedData(ex_policy)), terms_(encloser.options_.terms),
      omega_(encloser.options_.omega), pair_weight_regularization_(encloser.pair_weight_regularization_),
      reference_smoothing_length_(encloser.reference_smoothing_length_)
{
}

template <typename... Parameters>
inline void AphiApplyContactBlockDiagonalCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    const Real nu_i = nu_[index_i];
    const Vecd a_re_i = in_a_real_[index_i];
    const Vecd a_im_i = in_a_imag_[index_i];
    const Real phi_re_i = in_phi_real_[index_i];
    const Real phi_im_i = in_phi_imag_[index_i];

    Vecd lhs_a_real = Vecd::Zero();
    Vecd lhs_a_imag = Vecd::Zero();
    Real lhs_phi_real = Real(0);
    Real lhs_phi_imag = Real(0);

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
        const Real distance = r_ij_vec.norm();
        const Real distance_sq = r_ij_vec.squaredNorm();
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];

        if (terms_.laplace_phi)
        {
            const Real sigma_ij = AphiHarmonicMean(sigma_i, contact_sigma_[index_j]);
            const Real laplace_weight =
                sigma_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                             reference_smoothing_length_);
            lhs_phi_real += laplace_weight * phi_re_i;
            lhs_phi_imag += laplace_weight * phi_im_i;
        }

        if (terms_.laplace_a)
        {
            const Real nu_ij = AphiHarmonicMean(nu_i, contact_nu_[index_j]);
            const Real laplace_weight =
                nu_ij * AphiPairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                                          reference_smoothing_length_);
            lhs_a_real += laplace_weight * a_re_i;
            lhs_a_imag += laplace_weight * a_im_i;
        }

        if (terms_.grad_phi_coupling || terms_.div_sigma_a_coupling)
        {
            const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
            if (terms_.grad_phi_coupling)
            {
                lhs_a_real += sigma_i * g_ij * phi_re_i;
                lhs_a_imag += sigma_i * g_ij * phi_im_i;
            }
            if (terms_.div_sigma_a_coupling)
            {
                const Real sigma_ij = AphiHarmonicMean(sigma_i, contact_sigma_[index_j]);
                lhs_phi_real += omega_ * sigma_ij * g_ij.dot(a_im_i);
                lhs_phi_imag -= omega_ * sigma_ij * g_ij.dot(a_re_i);
            }
        }
    }

    out_a_real_[index_i] += lhs_a_real;
    out_a_imag_[index_i] += lhs_a_imag;
    out_phi_real_[index_i] += lhs_phi_real;
    out_phi_imag_[index_i] += lhs_phi_imag;
}

template <class ExecutionPolicy>
inline AphiApplyContactBlockDiagonalDynamicsBundle<ExecutionPolicy>::AphiApplyContactBlockDiagonalDynamicsBundle(
    SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation, const AphiBlockNames &input_block,
    const AphiBlockNames &output_block, const AphiMaterialNames &material_names, Real omega,
    const AphiLhsAssemblyOptions &options, Real pair_weight_regularization)
    : zero_output_(sph_body, output_block),
      apply_inner_(DynamicsArgs(inner_relation, input_block, output_block, material_names, omega, options,
                              pair_weight_regularization)),
      apply_contact_block_diagonal_(DynamicsArgs(contact_relation, input_block, output_block, material_names, omega,
                                                 options, pair_weight_regularization))
{
}

template <class ExecutionPolicy>
inline void AphiApplyContactBlockDiagonalDynamicsBundle<ExecutionPolicy>::exec()
{
    zero_output_.exec();
    apply_inner_.exec();
    apply_contact_block_diagonal_.exec();
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_MATRIX_FREE_OPERATOR_CK_HPP
