#ifndef FLUID_TIME_STEP_CK_HPP
#define FLUID_TIME_STEP_CK_HPP

#include "fluid_time_step_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy>
AcousticTimeStepCK::ReduceKernel::ReduceKernel(
    const ExecutionPolicy &ex_policy, AcousticTimeStepCK &encloser)
    : eos_(encloser.fluid_),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      p_(encloser.dv_p_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      force_(encloser.dv_force_->DelegatedDataField(ex_policy)),
      force_prior_(encloser.dv_force_prior_->DelegatedDataField(ex_policy)),
      h_min_(encloser.h_min_) {}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepSetup::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepSetup &encloser)
    : Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
AdvectionStepClose::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, AdvectionStepClose &encloser)
    : pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      dpos_(encloser.dv_dpos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_CK_HPP
