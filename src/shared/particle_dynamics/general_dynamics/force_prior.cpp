#include "force_prior.h"

namespace SPH
{
//=================================================================================================//
ForcePrior::ForcePrior(BaseParticles *base_particles, const std::string &force_name)
    : force_prior_(base_particles->registerSharedVariable<Vecd>("ForcePrior")),
      current_force_(base_particles->registerSharedVariable<Vecd>(force_name)),
      previous_force_(base_particles->registerSharedVariable<Vecd>("Previous" + force_name))
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
} // namespace SPH
