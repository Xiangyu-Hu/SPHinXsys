#include "fluid_integration.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
FluidInitialCondition::
    FluidInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), FluidDataSimple(sph_body),
      pos_(particles_->pos_), vel_(particles_->vel_) {}
//=================================================================================================//
ContinuumVolumeUpdate::ContinuumVolumeUpdate(SPHBody& sph_body)
    :LocalDynamics(sph_body), FluidDataSimple(sph_body),
    Vol_(particles_->Vol_), mass_(particles_->mass_), rho_(particles_->rho_) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
