#include "external_force.h"

namespace SPH
{
//=================================================================================================//
Gravity::Gravity(Vecd gravity_vector, Vecd reference_position)
    : reference_acceleration_(gravity_vector),
      zero_potential_reference_(reference_position) {}
//=================================================================================================//
Vecd Gravity::InducedAcceleration(const Vecd &position, Real physical_time)
{
    return reference_acceleration_;
}
//=================================================================================================//
Real Gravity::getPotential(const Vecd &position)
{
    return reference_acceleration_.dot(zero_potential_reference_ - position);
}
//=================================================================================================//
StartupAcceleration::StartupAcceleration(Vecd target_velocity, Real target_time)
    : Gravity(target_velocity / target_time), target_time_(target_time) {}
//=================================================================================================//
Vecd StartupAcceleration::InducedAcceleration(const Vecd &position, Real physical_time)
{
    Real time_factor = physical_time / target_time_;
    Real acceleration_factor = 0.5 * Pi * time_factor * sin(Pi * time_factor);

    return time_factor < 1.0 ? Gravity::InducedAcceleration() : Vecd::Zero();
}
//=================================================================================================//
Vecd IncreaseToFullGravity::InducedAcceleration(const Vecd &position, Real physical_time)
{
    Real time_factor = physical_time / time_to_full_gravity_;
    Vecd full_gravity = Gravity::InducedAcceleration();
    return time_factor < 1.0 ? time_factor * full_gravity : full_gravity;
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//