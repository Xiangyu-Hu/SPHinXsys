#ifndef RELAXATION_RESIDUE_CK_HPP
#define RELAXATION_RESIDUE_CK_HPP

#include "relaxation_residual_ck.h"

namespace SPH
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
RelaxationResidualBase<BaseInteractionType>::RelaxationResidualBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_residual_(this->particles_->template registerStateVariable<Vecd>("KernelGradientIntegral")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
RelaxationResidualCK<Inner<KernelCorrectionType, Parameters...>>::
    RelaxationResidualCK(Inner<Parameters...> &inner_relation)
    : BaseInteraction(inner_relation), kernel_correction_(this->particles_) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RelaxationResidualCK<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      residual_(encloser.dv_residual_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void RelaxationResidualCK<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Vecd residual = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);

        residual -= (this->correction_(index_i) + this->correction_(index_j)) * dW_ijV_j * e_ij;
    }
    residual_[index_i] = residual;
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
RelaxationResidualCK<Contact<Boundary, KernelCorrectionType, Parameters...>>::
    RelaxationResidualCK(Contact<Parameters...> &contact_relation)
    : BaseInteraction(contact_relation), kernel_correction_(this->particles_){}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RelaxationResidualCK<Contact<Boundary, KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)),
      residual_(encloser.dv_residual_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void RelaxationResidualCK<Contact<Boundary, KernelCorrectionType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Vecd residual = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);

        residual -= 2.0 * this->correction_(index_i) * dW_ijV_j * e_ij;
    }
    residual_[index_i] += residual;
}
//=================================================================================================//
} // namespace SPH
#endif // RELAXATION_RESIDUE_CK_HPP
