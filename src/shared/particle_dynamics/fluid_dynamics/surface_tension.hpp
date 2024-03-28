#pragma once

#include "surface_tension.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
SurfaceStressForce<DataDelegationType>::SurfaceStressForce(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_),
      color_gradient_(*this->particles_->template getVariableByName<Vecd>("ColorGradient")),
      surface_tension_force_(*this->particles_->template registerDiscreteVariable<Vecd>("SurfaceTensionForce")),
      surface_tension_stress_(*this->particles_->template getVariableByName<Matd>("SurfaceTensionStress")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
