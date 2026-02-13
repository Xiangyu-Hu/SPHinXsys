#pragma once
#include "density_correciton.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
DensitySummationPressure<Base, DataDelegationType>::DensitySummationPressure(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      rho_sum_(this->particles_->template registerStateVariableData<Real>("DensitySummation")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      rho0_(this->sph_body_->getBaseMaterial().ReferenceDensity()),
      inv_sigma0_(1.0 / this->getSPHAdaptation().LatticeNumberDensity()),
      W0_(this->getSPHAdaptation().getKernel()->W0(ZeroVecd)) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
