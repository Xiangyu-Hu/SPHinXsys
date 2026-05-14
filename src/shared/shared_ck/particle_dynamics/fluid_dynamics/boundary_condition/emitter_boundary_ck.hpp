#ifndef EMITTER_BOUNDARY_CK_HPP
#define EMITTER_BOUNDARY_CK_HPP

#include "emitter_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class OrientedBoxPartType, class ConditionFunction>
template <typename... Args>
EmitterInflowConditionCK<OrientedBoxPartType, ConditionFunction>::EmitterInflowConditionCK(
    OrientedBoxPartType &oriented_box_part, Args &&...args)
    : BaseLocalDynamics<OrientedBoxPartType>(oriented_box_part),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      inflow_velocity_(std::forward<Args>(args)...),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      sv_physical_time_(this->sph_system_->template getSystemVariableByName<Real>("PhysicalTime")) {}
//=================================================================================================//
template <class OrientedBoxPartType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
EmitterInflowConditionCK<OrientedBoxPartType, ConditionFunction>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      inflow_velocity_(encloser.inflow_velocity_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class OrientedBoxPartType, class ConditionFunction>
void EmitterInflowConditionCK<OrientedBoxPartType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    int aligned_axis = oriented_box_->ReferenceAxis();
    Transform &transform = oriented_box_->getTransform();
    Vecd frame_position = transform.shiftBaseStationToFrame(pos_[index_i]);
    Real current_axis_velocity = transform.xformBaseVecToFrame(vel_[index_i])[aligned_axis];
    Vecd frame_velocity = Vecd::Zero();
    frame_velocity[aligned_axis] =
        inflow_velocity_.getAxisVelocity(frame_position, current_axis_velocity, *physical_time_);
    vel_[index_i] = transform.xformFrameVecToBase(frame_velocity);
}
//=================================================================================================//
template <typename OrientedBoxPartType>
EmitterInflowInjectionCK<OrientedBoxPartType>::
    EmitterInflowInjectionCK(OrientedBoxPartType &oriented_box_part)
    : BaseLocalDynamics<OrientedBoxPartType>(oriented_box_part),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      spawn_real_particle_method_(this->particles_),
      rho0_(this->particles_->getBaseMaterial().ReferenceDensity()),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_p_(this->particles_->template getVariableByName<Real>("Pressure"))
{
    this->particles_->checkEnoughReserve();
}
//=================================================================================================//
template <typename OrientedBoxPartType>
EmitterInflowInjectionCK<OrientedBoxPartType>::FinishDynamics::
    FinishDynamics(EmitterInflowInjectionCK<OrientedBoxPartType> &encloser)
    : particles_(encloser.particles_) {}
//=================================================================================================//
template <typename OrientedBoxPartType>
void EmitterInflowInjectionCK<OrientedBoxPartType>::FinishDynamics::operator()()
{
    particles_->checkEnoughReserve();
}
//=================================================================================================//
template <typename OrientedBoxPartType>
template <class ExecutionPolicy, class EncloserType>
EmitterInflowInjectionCK<OrientedBoxPartType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      spawn_real_particle_(ex_policy, encloser.spawn_real_particle_method_),
      rho0_(encloser.rho0_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename OrientedBoxPartType>
void EmitterInflowInjectionCK<OrientedBoxPartType>::UpdateKernel::update(size_t index_i, Real dt)
{
    if (oriented_box_->checkUpperBound(pos_[index_i]))
    {
        Vecd original_position = pos_[index_i];
        spawn_real_particle_(index_i);
        pos_[index_i] = oriented_box_->getUpperPeriodic(original_position);
        rho_[index_i] = rho0_;
        p_[index_i] = 0.0;
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
WithinDisposerIndication::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      life_status_(encloser.dv_life_status_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void WithinDisposerIndication::UpdateKernel::update(size_t index_i, Real dt)
{
    if (oriented_box_->checkContain(pos_[index_i]) && index_i < *total_real_particles_)
    {
        life_status_[index_i] = 1; // mark as to delete but will not delete immediately
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // EMITTER_BOUNDARY_CK_HPP
