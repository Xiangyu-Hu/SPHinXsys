#include "base_material.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void Fluid::initializeLocalParameters(BaseParticles *base_particles)
{
    BaseMaterial::initializeLocalParameters(base_particles);
    base_particles->registerVariable(p_, "Pressure", getPressure(rho0_));
}
//=================================================================================================//
} // namespace SPH
