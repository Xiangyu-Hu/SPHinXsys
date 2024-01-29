#pragma once

#include "viscous_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
ViscousForce<DataDelegationType>::ViscousForce(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_), vel_(this->particles_->vel_),
      viscous_force_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForce")),
      mu_(DynamicCast<Fluid>(this, this->particles_->getBaseMaterial()).ReferenceViscosity()),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
