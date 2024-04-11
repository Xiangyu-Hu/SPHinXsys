#include "fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
FluidInitialCondition::
    FluidInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), FluidDataSimple(sph_body),
      pos_(particles_->ParticlePositions()),
      vel_(*particles_->getVariableByName<Vecd>("Velocity")) {}
//=================================================================================================//
ContinuumVolumeUpdate::ContinuumVolumeUpdate(SPHBody &sph_body)
    : LocalDynamics(sph_body), FluidDataSimple(sph_body),
      Vol_(particles_->VolumetricMeasures()),
      mass_(*particles_->getVariableByName<Real>("Mass")),
      rho_(*particles_->getVariableByName<Real>("Density")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
