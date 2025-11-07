#include "force_prior_ck.h"

namespace SPH
{
//=================================================================================================//
ForcePriorCK::ForcePriorCK(BaseParticles *particles, DiscreteVariable<Vecd> *dv_current_force)
    : dv_force_prior_(particles->registerStateVariable<Vecd>("ForcePrior")),
      dv_current_force_(dv_current_force),
      dv_previous_force_(particles->registerStateVariable<Vecd>("Previous" + dv_current_force->Name()))
{
    particles->addEvolvingVariable<Vecd>(dv_force_prior_);
    particles->addEvolvingVariable<Vecd>("Previous" + dv_current_force_->Name());
}
//=================================================================================================//
ForcePriorCK::ForcePriorCK(BaseParticles *particles, const std::string &force_name)
    : ForcePriorCK(particles, particles->template registerStateVariable<Vecd>(force_name)) {}
//=================================================================================================//
} // namespace SPH
