#ifndef APHI_LAPLACE_CK_HPP
#define APHI_LAPLACE_CK_HPP

#include "electromagnetic_dynamics/aphi_laplace_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
inline AphiPairwiseLaplaceCK<Base, DataType, RelationType<Parameters...>>::AphiPairwiseLaplaceCK(
    RelationType<Parameters...> &relation, const std::string &input_name,
    const std::string &coefficient_name, const std::string &output_name, Real pair_weight_regularization)
    : BaseInteraction(relation),
      input_name_(input_name),
      coefficient_name_(coefficient_name),
      output_name_(output_name),
      pair_weight_regularization_(pair_weight_regularization),
      reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      dv_input_(this->particles_->template getVariableByName<DataType>(input_name_)),
      dv_output_(this->particles_->template registerStateVariable<DataType>(output_name_, AphiZeroValue<DataType>::value())),
      dv_coefficient_(this->particles_->template getVariableByName<Real>(coefficient_name_))
{
    this->particles_->template addVariableToWrite<DataType>(output_name_);
}

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType, typename... Args>
inline AphiPairwiseLaplaceCK<Base, DataType, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
    : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      input_(encloser.dv_input_->DelegatedData(ex_policy)),
      output_(encloser.dv_output_->DelegatedData(ex_policy)),
      coefficient_(encloser.dv_coefficient_->DelegatedData(ex_policy)),
      pair_weight_regularization_(encloser.pair_weight_regularization_),
      reference_smoothing_length_(encloser.reference_smoothing_length_)
{
}

template <typename DataType, template <typename...> class RelationType, typename... Parameters>
inline void AphiPairwiseLaplaceCK<Base, DataType, RelationType<Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    (void)dt;
    const DataType value_i = input_[index_i];
    const Real coefficient_i = coefficient_[index_i];
    DataType laplace_i = AphiZeroValue<DataType>::value();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
        const Real distance = r_ij_vec.norm();
        const Real distance_sq = r_ij_vec.squaredNorm();
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Real harmonic_coefficient = AphiHarmonicMean(coefficient_i, coefficient_[index_j]);
        const Real pair_weight = harmonic_coefficient *
                                 AphiPairwiseNegativeLaplaceWeight(
                                     dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                                     reference_smoothing_length_);
        laplace_i += pair_weight * (value_i - input_[index_j]);
    }

    output_[index_i] = laplace_i;
}

template <typename DataType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiPairwiseLaplaceCK<Contact<DataType, Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser, contact_index),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_input_(encloser.dv_contact_input_[contact_index]->DelegatedData(ex_policy)),
      contact_coefficient_(encloser.dv_contact_coefficient_[contact_index]->DelegatedData(ex_policy))
{
}

template <typename DataType, typename... Parameters>
inline void AphiPairwiseLaplaceCK<Contact<DataType, Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const DataType value_i = this->input_[index_i];
    const Real coefficient_i = this->coefficient_[index_i];
    DataType laplace_i = AphiZeroValue<DataType>::value();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
        const Real distance = r_ij_vec.norm();
        const Real distance_sq = r_ij_vec.squaredNorm();
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Real harmonic_coefficient = AphiHarmonicMean(coefficient_i, contact_coefficient_[index_j]);
        const Real pair_weight = harmonic_coefficient *
                                 AphiPairwiseNegativeLaplaceWeight(
                                     dW_ijV_j, distance, distance_sq, this->pair_weight_regularization_,
                                     this->reference_smoothing_length_);
        laplace_i += pair_weight * (value_i - contact_input_[index_j]);
    }

    this->output_[index_i] += laplace_i;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_LAPLACE_CK_HPP
