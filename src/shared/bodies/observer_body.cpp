#include "observer_body.h"

#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
void ObserverBody::addObserverBodyToSPHSystem()
{
    sph_system_.addObservationBody(this);
}
//=================================================================================================//
} // namespace SPH
