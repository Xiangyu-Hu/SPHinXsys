#ifndef FORCE_PRIOR_HPP
#define FORCE_PRIOR_HPP

#include "force_prior.h"

namespace SPH
{
//=================================================================================================//
template <class GravityType>
GravityForce<GravityType>::GravityForce(SPHBody &sph_body, const GravityType &gravity)
    : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
      ForcePrior(particles_, "GravityForce"), gravity_(gravity),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      mass_(particles_->registerSharedVariable<Real>("Mass")),
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
} // namespace SPH
#endif // FORCE_PRIOR_HPP
