#ifndef FORCE_PRIOR_CK_HPP
#define FORCE_PRIOR_CK_HPP

#include "force_prior_ck.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
ForcePriorCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : force_prior_(encloser.dv_force_prior_->DelegatedData(ex_policy)),
      current_force_(encloser.dv_current_force_->DelegatedData(ex_policy)),
      previous_force_(encloser.dv_previous_force_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class GravityType>
GravityForceCK<GravityType>::GravityForceCK(SPHBody &sph_body, const GravityType &gravity)
    : LocalDynamics(sph_body), ForcePriorCK(this->particles_, "GravityForceCK"),
      gravity_(gravity),
      sv_physical_time_(sph_system_.getSystemVariableByName<Real>("PhysicalTime")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")) {}
//=================================================================================================//
template <class GravityType>
template <class ExecutionPolicy, class EncloserType>
GravityForceCK<GravityType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : ForcePriorCK::UpdateKernel(ex_policy, encloser),
      gravity_(encloser.gravity_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)) {}
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
