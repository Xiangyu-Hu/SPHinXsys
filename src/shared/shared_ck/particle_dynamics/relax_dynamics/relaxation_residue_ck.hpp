#ifndef RELAXATION_RESIDUE_CK_HPP
#define RELAXATION_RESIDUE_CK_HPP

#include "relaxation_residue_ck.h"

namespace SPH
{
//=================================================================================================//
template <class BaseInteractionType>
template <class DynamicsIdentifier>
RelaxationResidueBase<BaseInteractionType>::RelaxationResidueBase(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_residue_(this->particles_->template registerStateVariable<Vecd>("ZeroGradientResidue")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
RelaxationResidueCK<Inner<KernelCorrectionType, Parameters...>>::
    RelaxationResidueCK(Inner<Parameters...> &inner_relation)
    : BaseInteraction(inner_relation), kernel_correction_(this->particles_) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RelaxationResidueCK<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      correction_(ex_policy, encloser.kernel_correction_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      residue_(encloser.dv_residue_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename... Parameters>
void RelaxationResidueCK<Inner<KernelCorrectionType, Parameters...>>::InteractKernel::
    interact(size_t index_i, Real dt)
{
    Vecd residue = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
        const Vecd e_ij = this->e_ij(index_i, index_j);

        residue -= (this->correction_(index_i) + this->correction_(index_j)) * dW_ijV_j * e_ij;
    }
    residue_[index_i] = residue;
}
//=================================================================================================//
} // namespace SPH
#endif // RELAXATION_RESIDUE_CK_HPP
