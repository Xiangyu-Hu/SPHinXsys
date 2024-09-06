#ifndef FLUID_TIME_STEP_CK_HPP
#define FLUID_TIME_STEP_CK_HPP

#include "fluid_time_step_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class FluidType>
AcousticTimeStepCK<FluidType>::AcousticTimeStepCK(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      fluid_(DynamicCast<FluidType>(this, particles_->getBaseMaterial())),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure")),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")),
      h_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
template <class FluidType>
template <class ExecutionPolicy>
AcousticTimeStepCK<FluidType>::ReduceKernel::ReduceKernel(
    const ExecutionPolicy &ex_policy, AcousticTimeStepCK<FluidType> &encloser)
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedDataField(ex_policy)),
      h_min_(encloser.h_min_) {}
//=================================================================================================//
template <class FluidType>
Real AcousticTimeStepCK<FluidType>::ReduceKernel::reduce(size_t index_i, Real dt)
{
    Real acceleration_scale = 4.0 * h_min_ *
                              (force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i];
    return SMAX(eos_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm(), acceleration_scale);
}
//=================================================================================================//
template <class FluidType>
Real AcousticTimeStepCK<FluidType>::outputResult(Real reduced_value)
{
    // since the particle does not change its configuration in the acoustic time steps
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepSetup::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepSetup &encloser)
    : dv_Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      dv_mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      dv_dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepClose::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepClose &encloser)
    : dv_pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      dv_dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_CK_HPP
