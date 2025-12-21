#ifndef SURFACE_CORRECTION_HPP
#define SURFACE_CORRECTION_HPP

#include "surface_correction.h"

#include "level_set.hpp"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
LevelsetBounding::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      signed_distance_(ex_policy, encloser.level_set_, "LevelSet"),
      level_set_gradient_(ex_policy, encloser.level_set_, "LevelSetGradient"),
      constrained_distance_(encloser.constrained_distance_) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
LevelsetKernelGradientIntegral::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      residual_(encloser.dv_residual_->DelegatedData(ex_policy)),
      kernel_gradient_integral_(ex_policy, encloser.level_set_, "KernelGradient") {}
//=================================================================================================//
} // namespace SPH
#endif // SURFACE_CORRECTION_HPP