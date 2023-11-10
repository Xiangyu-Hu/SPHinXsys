#include "time_step_size_moving_velocity.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
AcousticTimeStepSizeMovingVelocity::AcousticTimeStepSizeMovingVelocity(SPHBody &sph_body, Real moving_velocity, Real acousticCFL)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      FluidDataSimple(sph_body), fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()), moving_velocity_(moving_velocity),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
Real AcousticTimeStepSizeMovingVelocity::reduce(size_t index_i, Real dt)
{
    return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + SMAX(vel_[index_i].norm(), moving_velocity_);
}
//=================================================================================================//
Real AcousticTimeStepSizeMovingVelocity::outputResult(Real reduced_value)
{
    // since the particle does not change its configuration in pressure relaxation step
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * smoothing_length_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
AdvectionTimeStepSizeMovingVelocity::AdvectionTimeStepSizeMovingVelocity(SPHBody &sph_body, Real U_ref, Real moving_velocity, Real advectionCFL)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_ref * U_ref),
      FluidDataSimple(sph_body), vel_(particles_->vel_),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL), moving_velocity_(moving_velocity),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()))
{
    Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
Real AdvectionTimeStepSizeMovingVelocity::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].squaredNorm();
}
//=================================================================================================//
Real AdvectionTimeStepSizeMovingVelocity::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * smoothing_length_min_ / (SMAX(moving_velocity_, speed_max, speed_ref_) + TinyReal);
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=====================================================================================================//
} // namespace SPH
  //=========================================================================================================//