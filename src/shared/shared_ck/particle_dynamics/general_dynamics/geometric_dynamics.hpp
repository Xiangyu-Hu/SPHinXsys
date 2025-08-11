#ifndef GEOMETRIC_DYNAMICS_HPP
#define GEOMETRIC_DYNAMICS_HPP

#include "execution_policy.h"
#include "geometric_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
InitialNearSurfaceInfo::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, InitialNearSurfaceInfo &encloser)
    : HostKernel(ex_policy, encloser),
      initial_shape_(encloser.initial_shape_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      n_(encloser.dv_n_->DelegatedData(ex_policy)),
      n0_(encloser.dv_n0_->DelegatedData(ex_policy)),
      phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
      phi0_(encloser.dv_phi0_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
SurfaceIndicationFromBodyShape::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, SurfaceIndicationFromBodyShape &encloser)
    : HostKernel(ex_policy, encloser),
      initial_shape_(encloser.initial_shape_),
      spacing_ref_(encloser.spacing_ref_),
      indicator_(encloser.dv_indicator_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
RandomizeParticlePositionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : HostKernel(ex_policy, encloser),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      randomize_scale_(encloser.randomize_scale_) {}
//=================================================================================================//
} // namespace SPH
#endif // GEOMETRIC_DYNAMICS_HPP
