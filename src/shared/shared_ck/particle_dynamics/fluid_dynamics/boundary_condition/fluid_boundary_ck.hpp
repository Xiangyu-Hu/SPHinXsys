#ifndef FLUID_BOUNDARY_CK_HPP
#define FLUID_BOUNDARY_CK_HPP

#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
InflowConditionCK<AlignedBoxPartType, ConditionFunction>::
    InflowConditionCK(AlignedBoxPartType &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      condition_function_(this->particles_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
InflowConditionCK<AlignedBoxPartType, ConditionFunction>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      condition_(ex_policy, encloser.condition_function_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
void InflowConditionCK<AlignedBoxPartType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    condition_(aligned_box_, index_i);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
EmitterInflowInjectionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      create_real_particle_(ex_policy, encloser.create_real_particle_method_),
      rho0_(encloser.rho0_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)) {}
//=================================================================================================//
void EmitterInflowInjectionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkUpperBound(pos_[index_i]))
    {
        create_real_particle_(index_i);
        pos_[index_i] = aligned_box_->getUpperPeriodic(pos_[index_i]); // Periodic bounding.
        rho_[index_i] = rho0_;
        p_[index_i] = 0.0;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_HPP
