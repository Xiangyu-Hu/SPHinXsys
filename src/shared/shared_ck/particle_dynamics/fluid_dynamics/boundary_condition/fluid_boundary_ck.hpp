#ifndef FLUID_BOUNDARY_CK_HPP
#define FLUID_BOUNDARY_CK_HPP

#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
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
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser),
    aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
    condition_(ex_policy, encloser.condition_function_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
void InflowConditionCK<AlignedBoxPartType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    condition_(aligned_box_, index_i);
}
//=================================================================================================//
EmitterInflowInjectionCK::
    EmitterInflowInjectionCK(AlignedBoxByParticle &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxByParticle>(aligned_box_part),
      buffer_(buffer), aligned_box_(aligned_box_part.svAlignedBox()),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure"))
{
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
EmitterInflowInjectionCK::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser),
    aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
    condition_(ex_policy, encloser.condition_function_) {}
//=================================================================================================//
void EmitterInflowInjectionCK::update(size_t original_index_i, Real dt)
{
    size_t sorted_index_i = sorted_id_[original_index_i];
    if (aligned_box_.checkUpperBound(pos_[sorted_index_i]))
    {
        mutex_switch_to_real_.lock();
        buffer_.checkEnoughBuffer(*particles_);
        particles_->createRealParticleFrom(sorted_index_i);
        mutex_switch_to_real_.unlock();

        /** Periodic bounding. */
        pos_[sorted_index_i] = aligned_box_.getUpperPeriodic(pos_[sorted_index_i]);
        rho_[sorted_index_i] = fluid_.ReferenceDensity();
        p_[sorted_index_i] = fluid_.getPressure(rho_[sorted_index_i]);
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_HPP
