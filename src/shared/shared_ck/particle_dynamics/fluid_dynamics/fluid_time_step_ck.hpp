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
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
template <class FluidType>
AcousticTimeStepCK<FluidType>::FinishDynamics::
    FinishDynamics(AcousticTimeStepCK<FluidType> &encloser)
    : h_min_(encloser.h_min_), acousticCFL_(encloser.acousticCFL_) {}
//=================================================================================================//
template <class FluidType>
Real AcousticTimeStepCK<FluidType>::FinishDynamics::Result(Real reduced_value)
{
    // since the particle does not change its configuration in the acoustic time steps
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
template <class FluidType>
template <class ExecutionPolicy>
AcousticTimeStepCK<FluidType>::ReduceKernel::ReduceKernel(
    const ExecutionPolicy &ex_policy, AcousticTimeStepCK<FluidType> &encloser)
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      h_min_(encloser.h_min_) {}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepSetup::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepSetup &encloser)
    : Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
UpdateParticlePosition::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, UpdateParticlePosition &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_CK_HPP
