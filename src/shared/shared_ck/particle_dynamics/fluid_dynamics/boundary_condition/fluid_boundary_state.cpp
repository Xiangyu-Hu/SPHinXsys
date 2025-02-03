#include "fluid_boundary_state.h"

namespace SPH
{
//=================================================================================================//
BaseStateCondition::BaseStateCondition(BaseParticles *particles)
    : dv_pos_(particles->getVariableByName<Vecd>("Position")),
      dv_vel_(particles->getVariableByName<Vecd>("Velocity")),
      dv_p_(particles->getVariableByName<Real>("Pressure")),
      dv_rho_(particles->getVariableByName<Real>("Density")) {}
//=================================================================================================//
} // namespace SPH

