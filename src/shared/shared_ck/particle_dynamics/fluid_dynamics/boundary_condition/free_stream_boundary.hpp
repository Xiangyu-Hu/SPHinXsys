#ifndef FLUID_BOUNDARY_STATE_HPP
#define FLUID_BOUNDARY_STATE_HPP

#include "fluid_boundary_state.h"

namespace SPH
{
//=================================================================================================//
template <typename ConditionFunction>
FreeStreamCondition<ConditionFunction>::FreeStreamCondition(
    SPHBody &sph_body, const ConditionFunction &free_stream_velocity)
    : LocalDynamics(sph_body), free_stream_velocity_(free_stream_velocity),
      rho0_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceDensity()),
      dv_rho_sum_(particles_->getVariableByName<Real>("DensitySummation")),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_indicator_(particles_->getVariableByName<int>("Indicator")),
      target_velocity(particles_),
      physical_time_(sph_system_->getSystemVariableByName<Real>("PhysicalTime")){};
//=================================================================================================//
template <typename ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
FreeStreamCondition<ConditionFunction>::ComputingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : free_stream_velocity_(encloser.free_stream_velocity_), rho0_(encloser.rho0_),
      rho_sum_(encloser.dv_rho_sum_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      physical_time_(encloser.physical_time_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename ConditionFunction>
void FreeStreamCondition<ConditionFunction>::ComputingKernel::update(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1)
    {
        Vecd current_velocity = vel_[index_i][0];
        Real target_velocity = free_stream_velocity_.getAxisVelocity(pos_[index_i], current_velocity, *physical_time_);
        Real corrected_velocity = current_velocity +
                                  (current_velocity - target_velocity) * SMIN(rho_sum_[index_i], rho0_) / rho0_;
        vel_[index_i][0] = corrected_velocity;
    }
};
//=================================================================================================//
} // namespace SPH
#endif // FLUID_BOUNDARY_STATE_HPP
