#ifndef KERNEL_INTEGRAL_HPP
#define KERNEL_INTEGRAL_HPP

#include "kernel_integral.h"

namespace SPH
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
KernelIntegralBase<BaseInteractionType>::KernelIntegralBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_kernel_integral_(
          this->particles_->template registerStateVariable<Real>("KernelIntegral")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
KernelIntegral<Inner<KernelCorrectionType, Parameters...>>::KernelIntegral(
    Inner<Parameters...> &inner_relation)
    : KernelIntegralBase<Interaction<Inner<Parameters...>>>(inner_relation) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
KernelIntegral<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      W0_(this->W(ZeroData<Vecd>::value)),
      kernel_integral_(encloser.dv_kernel_integral_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void KernelIntegral<Inner<KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real integral = W0_ * Vol_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        integral += this->W_ij(index_i, index_j) * Vol_[index_j];
    }
    kernel_integral_[index_i] = integral;
}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
KernelIntegral<Contact<KernelCorrectionType, Parameters...>>::
    KernelIntegral(Contact<Parameters...> &contact_relation)
    : KernelIntegralBase<Interaction<Contact<Parameters...>>>(contact_relation) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
KernelIntegral<Contact<KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      kernel_integral_(encloser.dv_kernel_integral_->DelegatedData(ex_policy)),
      contact_Vol_(encloser.dv_contact_Vol_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void KernelIntegral<Contact<KernelCorrectionType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real integral = 0.0;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        integral += this->W_ij(index_i, index_j) * contact_Vol_[index_j];
    }
    kernel_integral_[index_i] += integral;
}
//=================================================================================================//
} // namespace SPH
#endif // KERNEL_INTEGRAL_HPP