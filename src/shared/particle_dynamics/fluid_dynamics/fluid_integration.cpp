#include "fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
FluidInitialCondition::
    FluidInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), FluidDataSimple(sph_body),
      pos_(*base_particles_.getVariableByName<Vecd>("Position")),
      vel_(*particles_->registerSharedVariable<Vecd>("Velocity")) {}
//=================================================================================================//
ContinuumVolumeUpdate::ContinuumVolumeUpdate(SPHBody &sph_body)
    : LocalDynamics(sph_body), FluidDataSimple(sph_body),
      Vol_(*particles_->getVariableByName<Real>("VolumetricMeasure")),
      mass_(*particles_->getVariableByName<Real>("Mass")),
      rho_(*particles_->getVariableByName<Real>("Density")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
