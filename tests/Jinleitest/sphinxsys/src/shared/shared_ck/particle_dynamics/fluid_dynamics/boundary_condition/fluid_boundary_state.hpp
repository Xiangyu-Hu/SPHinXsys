#ifndef FLUID_BOUNDARY_STATE_HPP
#define FLUID_BOUNDARY_STATE_HPP

#include "fluid_boundary_state.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BaseStateCondition::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // FLUID_BOUNDARY_STATE_HPP
