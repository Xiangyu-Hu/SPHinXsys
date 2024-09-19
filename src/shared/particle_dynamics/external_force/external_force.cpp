#include "external_force.h"

namespace SPH
{
//=================================================================================================//
Gravity::Gravity(Vecd gravity_vector, Vecd reference_position)
    : reference_acceleration_(gravity_vector),
      zero_potential_reference_(reference_position) {}
//=================================================================================================//
StartupAcceleration::StartupAcceleration(Vecd target_velocity, Real target_time)
    : Gravity(target_velocity / target_time), target_time_(target_time) {}
//=================================================================================================//
IncreaseToFullGravity::IncreaseToFullGravity(Vecd gravity_vector, Real time_to_full_gravity)
    : Gravity(gravity_vector), time_to_full_gravity_(time_to_full_gravity) {}
//=================================================================================================//
} // namespace SPH
