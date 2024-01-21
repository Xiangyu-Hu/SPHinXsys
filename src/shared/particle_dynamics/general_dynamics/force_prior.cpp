#include "force_prior.h"

namespace SPH
{
//=================================================================================================//
ForcePrior::ForcePrior(BaseParticles *base_particles, const std::string &force_name)
    : force_prior_(base_particles->force_prior_),
      current_force_(*base_particles->registerSharedVariable<Vecd>(force_name)),
      previous_force_(*base_particles->registerSharedVariable<Vecd>("Previous" + force_name))
{
    base_particles->addVariableToRestart<Vecd>("Previous" + force_name);
}
//=================================================================================================//
void ForcePrior::update(size_t index_i, Real dt)
{
    force_prior_[index_i] += current_force_[index_i] - previous_force_[index_i];
    previous_force_[index_i] = current_force_[index_i];
}
//=================================================================================================//
GravityForce::GravityForce(SPHBody &sph_body, Gravity &gravity)
    : LocalDynamics(sph_body), ForcePrior(&base_particles_, "GravityForce"),
      GeneralDataDelegateSimple(sph_body), gravity_(gravity),
      pos_(base_particles_.pos_), mass_(base_particles_.mass_) {}
//=================================================================================================//
void GravityForce::update(size_t index_i, Real dt)
{
    current_force_[index_i] = mass_[index_i] * gravity_.InducedAcceleration(pos_[index_i]);
    ForcePrior::update(index_i, dt);
}
//=================================================================================================//
} // namespace SPH
