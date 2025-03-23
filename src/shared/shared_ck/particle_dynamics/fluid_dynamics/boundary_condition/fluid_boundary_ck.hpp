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
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      condition_function_(this->particles_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
InflowConditionCK<AlignedBoxPartType, ConditionFunction>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      condition_(ex_policy, encloser.condition_function_) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
void InflowConditionCK<AlignedBoxPartType, ConditionFunction>::UpdateKernel::update(size_t index_i, Real dt)
{
    condition_(aligned_box_, index_i, *physical_time_);
}
//=================================================================================================//
template <typename AlignedBoxPartType>
EmitterInflowInjectionCK<AlignedBoxPartType>::
    EmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      buffer_(buffer), sv_aligned_box_(aligned_box_part.svAlignedBox()),
      create_real_particle_method_(this->particles_),
      rho0_(this->particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_p_(this->particles_->template getVariableByName<Real>("Pressure"))
{
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
template <typename AlignedBoxPartType>
template <class ExecutionPolicy, class EncloserType>
EmitterInflowInjectionCK<AlignedBoxPartType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      create_real_particle_(ex_policy, encloser.create_real_particle_method_),
      rho0_(encloser.rho0_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)){}
//=================================================================================================//
template <typename AlignedBoxPartType>
void EmitterInflowInjectionCK<AlignedBoxPartType>::UpdateKernel::update(size_t index_i, Real dt)
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
template <typename AlignedBoxPartType, class ConditionFunction>
BufferEmitterInflowInjectionCK<AlignedBoxPartType, ConditionFunction>::
    BufferEmitterInflowInjectionCK(AlignedBoxPartType &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      buffer_(buffer),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      create_real_particle_method_(this->particles_),
      rho0_(this->particles_->getBaseMaterial().ReferenceDensity()),
      sound_speed_(0.0),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_p_(this->particles_->template getVariableByName<Real>("Pressure")),
      dv_buffer_particle_indicator_(this->particles_->template getVariableByName<int>("BufferParticleIndicator")),
      condition_function_(this->particles_),
      dv_previous_surface_indicator_(this->particles_->template getVariableByName<int>("PreviousSurfaceIndicator")),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      upper_bound_fringe_(0.5 * this->sph_body_.getSPHBodyResolutionRef())

{
    buffer_.checkParticlesReserved();
    Fluid &fluid_ = DynamicCast<Fluid>(this, this->particles_->getBaseMaterial());
    sound_speed_ = fluid_.getSoundSpeed(rho0_);
}
//=================================================================================================//
template <typename AlignedBoxPartType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
BufferEmitterInflowInjectionCK<AlignedBoxPartType, ConditionFunction>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(nullptr),
      create_real_particle_(ex_policy, encloser.create_real_particle_method_),
      rho0_(encloser.rho0_),
      sound_speed_(encloser.sound_speed_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)),
      condition_(ex_policy, encloser.condition_function_),
      previous_surface_indicator_(encloser.dv_previous_surface_indicator_->DelegatedData(ex_policy)),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      upper_bound_fringe_(encloser.upper_bound_fringe_)
{
    aligned_box_ = encloser.sv_aligned_box_->DelegatedData(ex_policy);
}
//=================================================================================================//
template <typename AlignedBoxPartType, class ConditionFunction>
void BufferEmitterInflowInjectionCK<AlignedBoxPartType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (!aligned_box_->checkInBounds(pos_[index_i]))
    {
        if (aligned_box_->checkUpperBound(pos_[index_i]))
        {
            if (buffer_particle_indicator_[index_i] == 1)
            {
                // if (index_i < this->particles_->TotalRealParticles())
                {
                    Vecd original_position = pos_[index_i];
                    size_t new_particle_index = create_real_particle_(index_i);
                    buffer_particle_indicator_[new_particle_index] = 0;
                    pos_[index_i] = aligned_box_->getUpperPeriodic(original_position);
                    p_[index_i] = condition_(index_i, *physical_time_);
                    rho_[index_i] = p_[index_i] / pow(sound_speed_, 2.0) + rho0_;
                    previous_surface_indicator_[index_i] = 1;
                }
            }
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
DisposerOutflowDeletionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      remove_real_particle_(ex_policy, encloser.remove_real_particle_method_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      is_deletable_(aligned_box_, pos_) {}
//=================================================================================================//
void DisposerOutflowDeletionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (is_deletable_(index_i))
    {
        remove_real_particle_(index_i, is_deletable_);
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TagBufferParticlesCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)) {}
//=================================================================================================//
void TagBufferParticlesCK::UpdateKernel::update(size_t index_i, Real dt)
{
    int buffer_indicator = 0;
    if (aligned_box_->checkInBounds(pos_[index_i]))
    {
        buffer_indicator = 1;
    }
    buffer_particle_indicator_[index_i] = buffer_indicator;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH

#endif // FLUID_BOUNDARY_CK_HPP
