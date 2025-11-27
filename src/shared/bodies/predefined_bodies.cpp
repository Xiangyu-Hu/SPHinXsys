#include "predefined_bodies.h"

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
void ObserverBody::addObserverBodyToSPHSystem()
{
    sph_system_.addObservationBody(this);
}
//=================================================================================================//
} // namespace SPH
