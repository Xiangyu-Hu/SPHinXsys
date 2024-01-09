#include "external_force.h"

namespace SPH
{
//=================================================================================================//
ExternalForce::ExternalForce() {}
//=================================================================================================//
Gravity::Gravity(Vecd global_acceleration, Vecd reference_position)
    : global_acceleration_(global_acceleration),
      zero_potential_reference_(reference_position),
      global_acceleration_device_(hostToDeviceVecd(global_acceleration_)),
      zero_potential_reference_device_(hostToDeviceVecd(zero_potential_reference_)) {}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//
