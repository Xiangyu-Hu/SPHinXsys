#ifndef EMITTER_BOUNDARY_CK_HPP
#define EMITTER_BOUNDARY_CK_HPP

#include "emitter_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>::
    EmitterInflowConditionCK(AlignedBoxPartType &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      condition_function_(this->particles_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      condition_(ex_policy, encloser.condition_function_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
void EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    condition_(aligned_box_, index_i, *physical_time_);
}
//=================================================================================================//
template <typename AlignedBoxPartType>
EmitterInflowInjectionCK<AlignedBoxPartType>::
    EmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      buffer_(buffer), sv_aligned_box_(aligned_box_part.svAlignedBox()),
      spawn_real_particle_method_(this->particles_),
      rho0_(this->particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_p_(this->particles_->template getVariableByName<Real>("Pressure"))
{
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
template <typename AlignedBoxPartType>
EmitterInflowInjectionCK<AlignedBoxPartType>::FinishDynamics::
    FinishDynamics(EmitterInflowInjectionCK<AlignedBoxPartType> &encloser)
    : particles_(encloser.particles_), buffer_(encloser.buffer_) {}
//=================================================================================================//
template <typename AlignedBoxPartType>
void EmitterInflowInjectionCK<AlignedBoxPartType>::FinishDynamics::operator()()
{
    buffer_.checkEnoughBuffer(*particles_);
}
//=================================================================================================//
template <typename AlignedBoxPartType>
template <class ExecutionPolicy, class EncloserType>
EmitterInflowInjectionCK<AlignedBoxPartType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      spawn_real_particle_(ex_policy, encloser.spawn_real_particle_method_),
      rho0_(encloser.rho0_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename AlignedBoxPartType>
void EmitterInflowInjectionCK<AlignedBoxPartType>::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkUpperBound(pos_[index_i]))
    {
        Vecd original_position = pos_[index_i];
        spawn_real_particle_(index_i);
        pos_[index_i] = aligned_box_->getUpperPeriodic(original_position);
        rho_[index_i] = rho0_;
        p_[index_i] = 0.0;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EMITTER_BOUNDARY_CK_HPP
