
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
RepulsionFactor<Contact<Parameters...>>::
    RepulsionFactor(Contact<Parameters...> &contact_relation)
    : BaseInteractionType(contact_relation, "RepulsionFactor")
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0 = this->contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0);
        dv_contact_mass_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Mass"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
RepulsionFactor<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
    : BaseInteractionType::InteractKernel(ex_policy, encloser, contact_index),
      repulsion_factor_(encloser.dv_repulsion_factor_->DelegatedData(ex_policy)),
      contact_inv_rho0_(encloser.contact_inv_rho0_[contact_index]),
      contact_mass_(encloser.dv_contact_mass_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void RepulsionFactor<Contact<Parameters...>>::InteractKernel::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * contact_inv_rho0_ * contact_mass_[index_j];
    }
    repulsion_factor_[index_i] = sigma;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_FACTOR_HPP
