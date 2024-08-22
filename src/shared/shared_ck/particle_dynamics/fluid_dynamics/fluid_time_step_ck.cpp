#include "fluid_time_step_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
AdvectionTimeStepCK::
    AdvectionTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      h_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")) {}
//=================================================================================================//
Real AdvectionTimeStepCK::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * h_min_ / (SMAX(speed_max, speed_ref_) + TinyReal);
}
//=================================================================================================//
AdvectionViscousTimeStepCK::AdvectionViscousTimeStepCK(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : AdvectionTimeStepCK(sph_body, U_ref, advectionCFL),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()))
{
    Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / h_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
