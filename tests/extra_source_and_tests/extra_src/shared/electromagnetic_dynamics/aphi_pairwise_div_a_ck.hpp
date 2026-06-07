#ifndef APHI_PAIRWISE_DIV_A_CK_HPP
#define APHI_PAIRWISE_DIV_A_CK_HPP

#include "electromagnetic_dynamics/aphi_pairwise_div_a_ck.h"
#include "electromagnetic_dynamics/aphi_laplace_ck.h"

namespace SPH
{
namespace electromagnetics
{

template <typename... Parameters>
inline AphiPairwiseVectorDivergenceCK<Inner<Parameters...>>::AphiPairwiseVectorDivergenceCK(
    Inner<Parameters...> &inner_relation, const std::string &input_name, const std::string &output_name)
    : BaseInteraction(inner_relation),
      dv_input_(this->particles_->template getVariableByName<Vecd>(input_name)),
      dv_output_(this->particles_->template getVariableByName<Real>(output_name))
{
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiPairwiseVectorDivergenceCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), input_(encloser.dv_input_->DelegatedData(ex_policy)),
      output_(encloser.dv_output_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiPairwiseVectorDivergenceCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd value_i = input_[index_i];
    Real divergence = Real(0);

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        divergence += g_ij.dot(value_i - input_[index_j]);
    }

    output_[index_i] = divergence;
}

template <typename... Parameters>
inline AphiPairwiseVectorDivergenceCK<Contact<Parameters...>>::AphiPairwiseVectorDivergenceCK(
    Contact<Parameters...> &contact_relation, const std::string &input_name, const std::string &output_name)
    : BaseInteraction(contact_relation),
      dv_input_(this->particles_->template getVariableByName<Vecd>(input_name)),
      dv_output_(this->particles_->template getVariableByName<Real>(output_name))
{
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_input_.push_back(contact_particles->template getVariableByName<Vecd>(input_name));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiPairwiseVectorDivergenceCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), input_(encloser.dv_input_->DelegatedData(ex_policy)),
      output_(encloser.dv_output_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_input_(encloser.dv_contact_input_[contact_index]->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiPairwiseVectorDivergenceCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd value_i = input_[index_i];
    Real divergence = Real(0);

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        divergence += g_ij.dot(value_i - contact_input_[index_j]);
    }

    output_[index_i] += divergence;
}

template <typename... Parameters>
inline AphiPairwiseScalarGradientCK<Inner<Parameters...>>::AphiPairwiseScalarGradientCK(
    Inner<Parameters...> &inner_relation, const std::string &input_name, const std::string &output_name)
    : BaseInteraction(inner_relation),
      dv_input_(this->particles_->template getVariableByName<Real>(input_name)),
      dv_output_(this->particles_->template getVariableByName<Vecd>(output_name))
{
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiPairwiseScalarGradientCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), input_(encloser.dv_input_->DelegatedData(ex_policy)),
      output_(encloser.dv_output_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiPairwiseScalarGradientCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real value_i = input_[index_i];
    Vecd gradient = Vecd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        gradient += g_ij * (value_i - input_[index_j]);
    }

    output_[index_i] = gradient;
}

template <typename... Parameters>
inline AphiPairwiseScalarGradientCK<Contact<Parameters...>>::AphiPairwiseScalarGradientCK(
    Contact<Parameters...> &contact_relation, const std::string &input_name, const std::string &output_name)
    : BaseInteraction(contact_relation),
      dv_input_(this->particles_->template getVariableByName<Real>(input_name)),
      dv_output_(this->particles_->template getVariableByName<Vecd>(output_name))
{
    for (auto *contact_particles : this->contact_particles_)
    {
        dv_contact_input_.push_back(contact_particles->template getVariableByName<Real>(input_name));
    }
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiPairwiseScalarGradientCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), input_(encloser.dv_input_->DelegatedData(ex_policy)),
      output_(encloser.dv_output_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      contact_input_(encloser.dv_contact_input_[contact_index]->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiPairwiseScalarGradientCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Real value_i = input_[index_i];
    Vecd gradient = Vecd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        gradient += g_ij * (value_i - contact_input_[index_j]);
    }

    output_[index_i] += gradient;
}

template <typename... Parameters>
inline AphiPairwiseVectorGradientCK<Inner<Parameters...>>::AphiPairwiseVectorGradientCK(
    Inner<Parameters...> &inner_relation, const std::string &input_name, const std::string &output_name)
    : BaseInteraction(inner_relation),
      dv_input_(this->particles_->template getVariableByName<Vecd>(input_name)),
      dv_output_(this->particles_->template getVariableByName<Matd>(output_name))
{
}

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline AphiPairwiseVectorGradientCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)), input_(encloser.dv_input_->DelegatedData(ex_policy)),
      output_(encloser.dv_output_->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void AphiPairwiseVectorGradientCK<Inner<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    const Vecd value_i = input_[index_i];
    Matd gradient = Matd::Zero();

    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd g_ij = AphiPairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
        const Vecd delta = value_i - input_[index_j];
        for (int m = 0; m != 3; ++m)
        {
            for (int n = 0; n != 3; ++n)
            {
                gradient(m, n) += g_ij[n] * delta[m];
            }
        }
    }

    output_[index_i] = gradient;
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_PAIRWISE_DIV_A_CK_HPP
