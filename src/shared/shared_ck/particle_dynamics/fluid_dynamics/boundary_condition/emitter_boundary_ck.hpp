#ifndef EMITTER_BOUNDARY_CK_HPP
#define EMITTER_BOUNDARY_CK_HPP

#include "emitter_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
template <typename... Args>
EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>::EmitterInflowConditionCK(
    AlignedBoxPartType &aligned_box_part, Args &&...args)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      inflow_velocity_(std::forward<Args>(args)...),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      sv_physical_time_(this->sph_system_->template getSystemVariableByName<Real>("PhysicalTime")) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      inflow_velocity_(encloser.inflow_velocity_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction>
void EmitterInflowConditionCK<AlignedBoxPartType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    int aligned_axis = aligned_box_->AlignmentAxis();
    Transform &transform = aligned_box_->getTransform();
    Vecd frame_position = transform.shiftBaseStationToFrame(pos_[index_i]);
    Real current_axis_velocity = transform.xformBaseVecToFrame(vel_[index_i])[aligned_axis];
    Vecd frame_velocity = Vecd::Zero();
    frame_velocity[aligned_axis] =
        inflow_velocity_.getAxisVelocity(frame_position, current_axis_velocity, *physical_time_);
    vel_[index_i] = transform.xformFrameVecToBase(frame_velocity);
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
template <class ExecutionPolicy, class EncloserType>
WithinDisposerIndication::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      life_status_(encloser.dv_life_status_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void WithinDisposerIndication::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]) && index_i < *total_real_particles_)
    {
            life_status_[index_i] = 1; // mark as to delete but will not delete immediately
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EMITTER_BOUNDARY_CK_HPP
