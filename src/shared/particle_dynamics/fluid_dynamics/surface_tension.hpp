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
    : ForcePrior(base_relation.getSPHBody(), "SurfaceTensionForce"),
      DataDelegationType(base_relation),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      color_gradient_(this->particles_->template getVariableDataByName<Vecd>("ColorGradient")),
      surface_tension_force_(this->particles_->template registerStateVariable<Vecd>("SurfaceTensionForce")),
      surface_tension_stress_(this->particles_->template getVariableDataByName<Matd>("SurfaceTensionStress")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
