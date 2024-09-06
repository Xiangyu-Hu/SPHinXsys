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
BaseForcePriorCK<DynamicsIdentifier>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy,
                               BaseForcePriorCK<DynamicsIdentifier> &encloser)
    : force_prior_(encloser.dv_force_prior_->DelegatedDataField(ex_policy)),
      current_force_(encloser.dv_current_force_->DelegatedDataField(ex_policy)),
      previous_force_(encloser.dv_previous_force_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
void BaseForcePriorCK<DynamicsIdentifier>::UpdateKernel::update(size_t index_i, Real dt)
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
GravityForceCK<GravityType>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, GravityForceCK<GravityType> &encloser)
    : ForcePriorCK::UpdateKernel(ex_policy, encloser),
      gravity_(encloser.gravity_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class GravityType>
void GravityForceCK<GravityType>::UpdateKernel::update(size_t index_i, Real dt)
{
    this->current_force_[index_i] =
        mass_[index_i] * gravity_.InducedAcceleration(pos_[index_i], *physical_time_);
    ForcePriorCK::UpdateKernel::update(index_i, dt);
}
//=================================================================================================//
} // namespace SPH
#endif // FORCE_PRIOR_CK_HPP
