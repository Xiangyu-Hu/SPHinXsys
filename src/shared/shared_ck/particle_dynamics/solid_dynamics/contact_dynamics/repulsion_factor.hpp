
#ifndef REPULSION_FACTOR_HPP
#define REPULSION_FACTOR_HPP

#include "base_particles.hpp"
#include "repulsion_factor.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
RepulsionFactor<Base, Contact<Parameters...>>::RepulsionFactor(
    DynamicsIdentifier &identifier, const std::string &factor_name)
    : BaseInteractionType(identifier),
      dv_repulsion_factor_(this->particles_->template registerStateVariable<Real>(factor_name)) {}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
RepulsionFactor<Contact<Parameters...>>::
    RepulsionFactor(DynamicsIdentifier &identifier)
    : BaseInteractionType(identifier, "RepulsionFactor"),
    dv_contact_Vol_ref_(this->contact_particles_->template registerStateVariableFrom<Real>(
        "VolumetricMeasureRef", "VolumetricMeasure")){}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RepulsionFactor<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseInteractionType::InteractKernel(ex_policy, encloser),
      repulsion_factor_(encloser.dv_repulsion_factor_->DelegatedData(ex_policy)),
      contact_Vol_ref_(encloser.dv_contact_Vol_ref_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void RepulsionFactor<Contact<Parameters...>>::InteractKernel::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * contact_Vol_ref_[index_j];
    }
    repulsion_factor_[index_i] = sigma;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_FACTOR_HPP
