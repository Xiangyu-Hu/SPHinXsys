#ifndef FORCE_PRIOR_CK_HPP
#define FORCE_PRIOR_CK_HPP

#include "force_prior_ck.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier>
BaseForcePriorCK<DynamicsIdentifier>::
    BaseForcePriorCK(DynamicsIdentifier &identifier, const std::string &force_name)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_force_prior_(this->particles_->template registerStateVariableOnly<Vecd>("ForcePrior")),
      dv_current_force_(this->particles_->template registerStateVariableOnly<Vecd>(force_name)),
      dv_previous_force_(this->particles_->template registerStateVariableOnly<Vecd>("Previous" + force_name))
{
    this->particles_->template addVariableToRestart<Vecd>("Previous" + force_name);
    this->particles_->template addVariableToSort<Vecd>("Previous" + force_name);
}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy>
BaseForcePriorCK<DynamicsIdentifier>::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    BaseForcePriorCK<DynamicsIdentifier> &base_force_prior)
    : force_prior_(base_force_prior.dv_force_prior_->DelegatedDataField(ex_policy)),
      current_force_(base_force_prior.dv_current_force_->DelegatedDataField(ex_policy)),
      previous_force_(base_force_prior.dv_previous_force_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy>
void BaseForcePriorCK<DynamicsIdentifier>::
    ComputingKernel<ExecutionPolicy>::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += current_force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = current_force_[index_i];
}
//=================================================================================================//
template <class GravityType>
GravityForceCK<GravityType>::GravityForceCK(SPHBody &sph_body, const GravityType &gravity)
    : ForcePriorCK(sph_body, "GravityForceCK"), gravity_(gravity),
      sv_physical_time_(sph_system_.getSystemVariableByName<Real>("PhysicalTime")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")) {}
//=================================================================================================//
template <class GravityType>
template <class ExecutionPolicy>
GravityForceCK<GravityType>::ComputingKernel<ExecutionPolicy>::ComputingKernel(
    const ExecutionPolicy &ex_policy, GravityForceCK<GravityType> &gravity_force)
    : ForcePriorCK::ComputingKernel<ExecutionPolicy>(ex_policy, gravity_force),
      gravity_(gravity_force.gravity_),
      physical_time_(gravity_force.sv_physical_time_->DelegatedData(ex_policy)),
      pos_(gravity_force.dv_pos_->DelegatedDataField(ex_policy)),
      mass_(gravity_force.dv_mass_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class GravityType>
template <class ExecutionPolicy>
void GravityForceCK<GravityType>::ComputingKernel<ExecutionPolicy>::update(size_t index_i, Real dt)
{
    this->current_force_[index_i] =
        mass_[index_i] * gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
    ForcePriorCK::ComputingKernel<ExecutionPolicy>::update(index_i, dt);
}
//=================================================================================================//
} // namespace SPH
#endif // FORCE_PRIOR_CK_HPP
