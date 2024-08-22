#ifndef FLUID_TIME_STEP_CK_HPP
#define FLUID_TIME_STEP_CK_HPP

#include "fluid_time_step_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionTimeStepCK::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy, AdvectionTimeStepCK &advection_time_step)
    : h_min_(advection_time_step.h_min_),
      mass_(advection_time_step.dv_mass_->DelegatedDataField(ex_policy)),
      vel_(advection_time_step.dv_vel_->DelegatedDataField(ex_policy)),
      force_(advection_time_step.dv_force_->DelegatedDataField(ex_policy)),
      force_prior_(advection_time_step.dv_force_prior_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
Real AdvectionTimeStepCK::ComputingKernel<ExecutionPolicy>::reduce(size_t index_i, Real dt)
{
    Real acceleration_scale = 4.0 * h_min_ *
                              (force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i];
    return SMAX(vel_[index_i].squaredNorm(), acceleration_scale);
};
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionViscousTimeStepCK::ComputingKernel<ExecutionPolicy>::ComputingKernel(
    const ExecutionPolicy &ex_policy, AdvectionViscousTimeStepCK &advection_viscous_time_step)
    : AdvectionTimeStepCK::ComputingKernel<ExecutionPolicy>(ex_policy, advection_viscous_time_step){};
//=================================================================================================//
template <class ExecutionPolicy>
Real AdvectionViscousTimeStepCK::ComputingKernel<ExecutionPolicy>::reduce(size_t index_i, Real dt)
{
    return AdvectionTimeStepCK::ComputingKernel<ExecutionPolicy>::reduce(index_i, dt);
};
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_CK_HPP
