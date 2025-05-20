#ifndef FORCE_PRIOR_HPP
#define FORCE_PRIOR_HPP

#include "force_prior.h"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier>
BaseForcePrior<DynamicsIdentifier>::
    BaseForcePrior(DynamicsIdentifier &identifier, const std::string &force_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      force_prior_(this->particles_->template registerStateVariable<Vecd>("ForcePrior")),
      current_force_(this->particles_->template registerStateVariable<Vecd>(force_name)),
      previous_force_(this->particles_->template registerStateVariable<Vecd>("Previous" + force_name))
{
    this->particles_->template addEvolvingVariable<Vecd>("Previous" + force_name);
}
//=================================================================================================//
template <class DynamicsIdentifier>
void BaseForcePrior<DynamicsIdentifier>::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += current_force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = current_force_[index_i];
}
//=================================================================================================//
template <class GravityType>
GravityForce<GravityType>::GravityForce(SPHBody &sph_body, const GravityType &gravity)
    : ForcePrior(sph_body, "GravityForce"), gravity_(gravity),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      mass_(particles_->registerStateVariable<Real>("Mass")),
      physical_time_(sph_system_.getSystemVariableDataByName<Real>("PhysicalTime")) {}
//=================================================================================================//
template <class GravityType>
void GravityForce<GravityType>::GravityForce::update(size_t index_i, Real dt)
{
    current_force_[index_i] =
        mass_[index_i] * gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
    ForcePrior::update(index_i, dt);
}
//=================================================================================================//
} // namespace SPH
#endif // FORCE_PRIOR_HPP
