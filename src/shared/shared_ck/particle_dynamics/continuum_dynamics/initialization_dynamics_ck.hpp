#ifndef INITIALIZATION_DYNAMICS_CK_HPP
#define INITIALIZATION_DYNAMICS_CK_HPP

#include "initialization_dynamics_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
ContinuumInitialConditionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      stress_tensor_3D_(encloser.dv_stress_tensor_3D_->DelegatedData(ex_policy)){};
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
#endif // INITIALIZATION_DYNAMICS_CK_HPP