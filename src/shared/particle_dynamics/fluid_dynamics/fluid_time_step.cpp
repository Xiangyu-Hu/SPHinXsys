#include "fluid_time_step.h"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      FluidDataSimple(sph_body), fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
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
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_ref * U_ref),
      FluidDataSimple(sph_body), mass_(particles_->mass_), vel_(particles_->vel_),
      force_(particles_->force_), force_prior_(particles_->force_prior_),
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

SRDViscousTimeStepSize::SRDViscousTimeStepSize(SPHBody &sph_body, Real diffusionCFL) : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
                                                                                       FluidDataSimple(sph_body),
                                                                                       smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
                                                                                       rho_(this->particles_->rho_),
                                                                                       mu_srd_(*this->particles_->getVariableByName<Real>("SRDViscosity")),
                                                                                       diffusionCFL(diffusionCFL)
{
}

Real SRDViscousTimeStepSize::outputResult(size_t index_i, Real dt)
{
    return this->diffusionCFL * smoothing_length_ * smoothing_length_ * rho_ / max_viscosity;
}

Real SRDViscousTimeStepSize::reduce(size_t index_i, Real dt)
{
    max_viscosity = SMAX(mu_srd_[index_i], max_viscosity);
    return max_viscosity;
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=====================================================================================================//
} // namespace SPH
  //=========================================================================================================//