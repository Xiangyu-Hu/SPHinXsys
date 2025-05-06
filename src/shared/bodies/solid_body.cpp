#include "solid_body.h"

#include "base_material.h"
#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
void SolidBody::addSolidBodyToSPHSystem()
{
    sph_system_.addSolidBody(this);
}
//=================================================================================================//
} // namespace SPH
