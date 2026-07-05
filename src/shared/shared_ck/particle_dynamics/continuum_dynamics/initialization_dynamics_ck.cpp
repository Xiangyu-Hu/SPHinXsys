#ifndef INITIALIZATION_DYNAMICS_CK_HPP
#define INITIALIZATION_DYNAMICS_CK_HPP

#include "initialization_dynamics_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
  ContinuumInitialConditionCK::ContinuumInitialConditionCK(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_pos_(this->particles_->template registerStateVariable<Vecd>("Position")),
      dv_vel_(this->particles_->template registerStateVariable<Vecd>("Velocity")),
      dv_stress_tensor_3D_(this->particles_->template registerStateVariable<Mat3d>("StressTensor3D")){}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
#endif // INITIALIZATION_DYNAMICS_CK_HPP