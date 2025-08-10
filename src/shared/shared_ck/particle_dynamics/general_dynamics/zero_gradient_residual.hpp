#ifndef ZERO_GRADIENT_RESIDUAL_HPP
#define ZERO_GRADIENT_RESIDUAL_HPP

#include "zero_gradient_residual.h"

namespace SPH
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
ZeroGradientResidualBase<BaseInteractionType>::ZeroGradientResidualBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_zero_gradient_residue_(
          this->particles_->template registerStateVariable<Vecd>("ZeroGradientResidue")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
ZeroGradientResidual<Inner<KernelCorrectionType, Parameters...>>::ZeroGradientResidual(
    Inner<Parameters...> &inner_relation)
    : ZeroGradientResidualBase<Interaction<Inner<Parameters...>>>(inner_relation),
      kernel_correction_(this->particles_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrectionType must derive from KernelCorrection!");
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ZeroGradientResidual<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void ZeroGradientResidual<Inner<KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd inconsistency = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);
        inconsistency -= (correction_(index_i) + correction_(index_j)) * dW_ijV_j * e_ij;
    }
    zero_gradient_residue_[index_i] = inconsistency;
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
ZeroGradientResidual<Contact<Boundary, KernelCorrectionType, Parameters...>>::
    ZeroGradientResidual(Contact<Parameters...> &contact_relation)
    : ZeroGradientResidualBase<Interaction<Contact<Parameters...>>>(contact_relation),
      kernel_correction_(this->particles_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrectionType must derive from KernelCorrection!");
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
ZeroGradientResidual<Contact<Boundary, KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void ZeroGradientResidual<Contact<Boundary, KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd inconsistency = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * contact_Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);
        inconsistency -= 2.0 * correction_(index_i) * dW_ijV_j * e_ij;
    }
    zero_gradient_residue_[index_i] += inconsistency;
}
//=================================================================================================//
} // namespace SPH
#endif // ZERO_GRADIENT_RESIDUAL_HPP