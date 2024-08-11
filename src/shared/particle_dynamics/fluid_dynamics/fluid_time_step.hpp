#ifndef FLUID_TIME_STEP_HPP
#define FLUID_TIME_STEP_HPP

#include "fluid_time_step.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionTimeStep::AdvectionTimeStep(const ExecutionPolicy &exec_policy,
                                     SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      mass_(particles_->getVariableDataByName<Real>(exec_policy, "Mass")),
      vel_(particles_->getVariableDataByName<Vecd>(exec_policy, "Velocity")),
      force_(particles_->getVariableDataByName<Vecd>(exec_policy, "Force")),
      force_prior_(particles_->getVariableDataByName<Vecd>(exec_policy, "ForcePrior")),
      h_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL) {}
//=================================================================================================//
template <typename... Args>
AdvectionViscousTimeStep::AdvectionViscousTimeStep(Args &...args)
    : AdvectionTimeStep(std::forward<Args>(args)...),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()))
{
    Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / h_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_HPP
