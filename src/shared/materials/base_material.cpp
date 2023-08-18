#include "base_material.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void Fluid::registerReloadLocalParameters(BaseParticles *base_particles)
{
    base_particles->registerVariable(p_, "Pressure");
}
//=================================================================================================//
} // namespace SPH
