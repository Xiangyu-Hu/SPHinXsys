#ifndef FORCE_PRIOR_HPP
#define FORCE_PRIOR_HPP

#include "force_prior.h"

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
    this->particles_->template addVariableToRestart<Vecd>("Previous" + force_name);
    this->particles_->template addVariableToSort<Vecd>("Previous" + force_name);
}
//=================================================================================================//
template <class DynamicsIdentifier>
void BaseForcePrior<DynamicsIdentifier>::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += current_force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = current_force_[index_i];
}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy>
BaseForcePrior<DynamicsIdentifier>::
    BaseForcePrior(const ExecutionPolicy &execution_policy,
                   DynamicsIdentifier &identifier, const std::string &force_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      force_prior_(this->particles_->template registerStateVariable<Vecd>(execution_policy, "ForcePrior")),
      current_force_(this->particles_->template registerStateVariable<Vecd>(execution_policy, force_name)),
      previous_force_(this->particles_->template registerStateVariable<Vecd>(execution_policy, "Previous" + force_name))
{
    this->particles_->template addVariableToRestart<Vecd>("Previous" + force_name);
    this->particles_->template addVariableToSort<Vecd>("Previous" + force_name);
}
//=================================================================================================//
template <class DynamicsIdentifier>
BaseForcePrior<DynamicsIdentifier>::ComputingKernel::
    ComputingKernel(BaseForcePrior<DynamicsIdentifier> &base_force_prior)
    : force_prior_(base_force_prior.force_prior_),
      current_force_(base_force_prior.current_force_),
      previous_force_(base_force_prior.previous_force_) {}
//=================================================================================================//
template <class DynamicsIdentifier>
void BaseForcePrior<DynamicsIdentifier>::ComputingKernel::update(size_t index_i, Real dt)
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
    current_force_[index_i] = mass_[index_i] *
                              gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
    ForcePrior::update(index_i, dt);
}
//=================================================================================================//
template <class GravityType>
template <class ExecutionPolicy>
GravityForce<GravityType>::
    GravityForce(const ExecutionPolicy &execution_policy,
                 SPHBody &sph_body, const GravityType &gravity)
    : ForcePrior(execution_policy, sph_body, "GravityForce"), gravity_(gravity),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      mass_(particles_->getVariableDataByName<Real>(execution_policy, "Mass")),
      v_physical_time_(sph_system_.getSystemVariableByName<Real>(execution_policy, "PhysicalTime")) {}
//=================================================================================================//
template <class GravityType>
GravityForce<GravityType>::ComputingKernel::ComputingKernel(GravityForce<GravityType> &gravity_force)
    : ForcePrior::ComputingKernel(gravity_force),
      gravity_(gravity_force.gravity_), pos_(gravity_force.pos_), mass_(gravity_force.mass_),
      physical_time_(gravity_force.v_physical_time_->ValueAddress()) {}
//=================================================================================================//
template <class GravityType>
void GravityForce<GravityType>::ComputingKernel::update(size_t index_i, Real dt)
{
    current_force_[index_i] = mass_[index_i] *
                              gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
    ForcePrior::ComputingKernel::update(index_i, dt);
}
//=================================================================================================//
} // namespace SPH
#endif // FORCE_PRIOR_HPP
