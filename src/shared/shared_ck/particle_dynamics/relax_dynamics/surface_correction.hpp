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
      signed_distance_(encloser.level_set_.getProbeSignedDistance(ex_policy)),
      normal_direction_(encloser.level_set_.getProbeNormalDirection(ex_policy)),
      constrained_distance_(encloser.constrained_distance_) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
LevelsetKernelGradientIntegral::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      residual_(encloser.dv_residual_->DelegatedData(ex_policy)),
      kernel_gradient_integral_(encloser.level_set_.getProbeKernelGradientIntegral(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // SURFACE_CORRECTION_HPP