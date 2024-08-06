#include "force_prior.h"

namespace SPH
{
//=================================================================================================//
ForcePrior::ForcePrior(BaseParticles *base_particles, const std::string &force_name)
    : force_prior_(*base_particles->registerSharedVariable<Vecd>("ForcePrior")),
      current_force_(*base_particles->registerSharedVariable<Vecd>(force_name)),
      previous_force_(*base_particles->registerSharedVariable<Vecd>("Previous" + force_name))
{
    base_particles->addVariableToRestart<Vecd>("Previous" + force_name);
    base_particles->addVariableToSort<Vecd>("Previous" + force_name);
}
//=================================================================================================//
void ForcePrior::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += current_force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = current_force_[index_i];
}
//=================================================================================================//
GravityForce::GravityForce(SPHBody &sph_body, Gravity &gravity)
    : LocalDynamics(sph_body), DataDelegateSimple(sph_body),
      ForcePrior(particles_, "GravityForce"), gravity_(gravity),
      pos_(*particles_->getVariableDataByName<Vecd>("Position")),
      mass_(*particles_->registerSharedVariable<Real>("Mass")) {}
//=================================================================================================//
void GravityForce::update(size_t index_i, Real dt)
{
    current_force_[index_i] = mass_[index_i] * gravity_.InducedAcceleration(pos_[index_i]);
    ForcePrior::update(index_i, dt);
}
//=================================================================================================//
} // namespace SPH
