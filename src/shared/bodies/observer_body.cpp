#include "observer_body.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
void ObserverBody::addObserverBodyToSPHSystem()
{
    sph_system_.observation_bodies_.push_back(this);
}
//=================================================================================================//
} // namespace SPH
