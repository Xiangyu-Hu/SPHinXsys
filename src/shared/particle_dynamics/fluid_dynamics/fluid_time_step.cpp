#include "fluid_time_step.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      DataDelegateSimple(sph_body), fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      rho_(*particles_->getVariableDataByName<Real>("Density")),
      p_(*particles_->getVariableDataByName<Real>("Pressure")),
      mass_(*particles_->getVariableDataByName<Real>("Mass")),
      vel_(*particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(*particles_->getVariableDataByName<Vecd>("Force")),
      force_prior_(*particles_->getVariableDataByName<Vecd>("ForcePrior")),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    Real acceleration_scale = 4.0 * smoothing_length_min_ *
                              (force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i];
    return SMAX(fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm(), acceleration_scale);
}
//=================================================================================================//
Real AcousticTimeStepSize::outputResult(Real reduced_value)
{
    // since the particle does not change its configuration in pressure relaxation step
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * smoothing_length_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
AdvectionTimeStepSizeForImplicitViscosity::
    AdvectionTimeStepSizeForImplicitViscosity(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      DataDelegateSimple(sph_body),
      mass_(*particles_->getVariableDataByName<Real>("Mass")),
      vel_(*particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(*particles_->getVariableDataByName<Vecd>("Force")),
      force_prior_(*particles_->getVariableDataByName<Vecd>("ForcePrior")),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL) {}
//=================================================================================================//
Real AdvectionTimeStepSizeForImplicitViscosity::reduce(size_t index_i, Real dt)
{
    Real acceleration_scale = 4.0 * smoothing_length_min_ *
                              (force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i];
    return SMAX(vel_[index_i].squaredNorm(), acceleration_scale);
}
//=================================================================================================//
Real AdvectionTimeStepSizeForImplicitViscosity::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * smoothing_length_min_ / (SMAX(speed_max, speed_ref_) + TinyReal);
}
//=================================================================================================//
AdvectionTimeStepSize::AdvectionTimeStepSize(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : AdvectionTimeStepSizeForImplicitViscosity(sph_body, U_ref, advectionCFL),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()))
{
    Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
Real AdvectionTimeStepSize::reduce(size_t index_i, Real dt)
{
    return AdvectionTimeStepSizeForImplicitViscosity::reduce(index_i, dt);
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=====================================================================================================//
} // namespace SPH
  //=========================================================================================================//