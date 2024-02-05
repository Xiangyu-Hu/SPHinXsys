#pragma once

#include "non_newtonian_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class DataDelegationType>
template <class BaseRelationType>
HerschelBulkleyAcceleration<DataDelegationType>::HerschelBulkleyAcceleration(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_), vel_(this->particles_->vel_),
      viscous_force_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForce")),
      // combined_velocity_gradient_(this->particles_->getVariableByName<Matd>("CombinedVelocityGradient")),
      herschel_bulkley_fluid_(DynamicCast<HerschelBulkleyFluid>(this, this->particles_->getBaseMaterial())),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength())
{
}
} // namespace fluid_dynamics
} // namespace SPH