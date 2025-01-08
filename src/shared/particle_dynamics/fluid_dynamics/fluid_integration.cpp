#include "fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
FluidInitialCondition::
    FluidInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariable<Vecd>("Velocity")) {}
//=================================================================================================//
ContinuumVolumeUpdate::ContinuumVolumeUpdate(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      rho_(particles_->getVariableDataByName<Real>("Density")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
