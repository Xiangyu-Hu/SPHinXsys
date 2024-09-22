#ifndef GEOMETRIC_DYNAMICS_HPP
#define GEOMETRIC_DYNAMICS_HPP

#include "execution_policy.h"
#include "geometric_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
NormalFromBodyShapeCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, NormalFromBodyShapeCK &encloser)
    : shape_(encloser.shape_),
      pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)),
      n_(encloser.dv_n_->DelegatedDataField(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedDataField(ex_policy)),
      phi_(encloser.dv_phi_->DelegatedDataField(ex_policy)),
      phi0_(encloser.dv_phi0_->DelegatedDataField(ex_policy))
{
    // not implemented for device policy due to virtual function call in inital_shape_,
    // which is not allowed in device code
    static_assert(!std::is_base_of<execution::ParallelDevicePolicy, ExecutionPolicy>::value,
                  "This compute kernel is not designed for execution::ParallelDevicePolicy!");
}
//=================================================================================================//
} // namespace SPH
#endif // GEOMETRIC_DYNAMICS_HPP
