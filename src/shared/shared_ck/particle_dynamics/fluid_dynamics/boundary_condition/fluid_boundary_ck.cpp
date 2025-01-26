#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
BaseStateCondition::BaseStateCondition(BaseParticles *particles)
    : dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")),
      dv_rho_(particles_->getVariableByName<Real>("Density")) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_HPP
