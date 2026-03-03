#ifndef KERNEL_GRADIENT_INTEGRAL_HPP
#define KERNEL_GRADIENT_INTEGRAL_HPP

#include "kernel_gradient_integral.h"

namespace SPH
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
KernelGradientIntegralBase<BaseInteractionType>::KernelGradientIntegralBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_kernel_gradient_integral_(
          this->particles_->template registerStateVariable<Vecd>("KernelGradientIntegral")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
KernelGradientIntegral<Inner<KernelCorrectionType, Parameters...>>::KernelGradientIntegral(
    Inner<Parameters...> &inner_relation)
    : KernelGradientIntegralBase<Interaction<Inner<Parameters...>>>(inner_relation),
      kernel_correction_(this->particles_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrectionType must derive from KernelCorrection!");
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
KernelGradientIntegral<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      kernel_gradient_integral_(encloser.dv_kernel_gradient_integral_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void KernelGradientIntegral<Inner<KernelCorrectionType, Parameters...>>::
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
    kernel_gradient_integral_[index_i] = inconsistency;
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
KernelGradientIntegral<Contact<Boundary, KernelCorrectionType, Parameters...>>::
    KernelGradientIntegral(Contact<Parameters...> &contact_relation)
    : KernelGradientIntegralBase<Interaction<Contact<Parameters...>>>(contact_relation),
      kernel_correction_(this->particles_)
{
    static_assert(std::is_base_of<KernelCorrection, KernelCorrectionType>::value,
                  "KernelCorrectionType must derive from KernelCorrection!");
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
KernelGradientIntegral<Contact<Boundary, KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      correction_(ex_policy, encloser.kernel_correction_),
      kernel_gradient_integral_(encloser.dv_kernel_gradient_integral_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void KernelGradientIntegral<Contact<Boundary, KernelCorrectionType, Parameters...>>::
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
    kernel_gradient_integral_[index_i] += inconsistency;
}
//=================================================================================================//
} // namespace SPH
#endif // KERNEL_GRADIENT_RESIDUAL_HPP