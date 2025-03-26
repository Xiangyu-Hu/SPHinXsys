#ifndef PRESSURE_BOUNDARY_CK_HPP
#define PRESSURE_BOUNDARY_CK_HPP

#include "pressure_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TagBufferParticlesCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : part_id_(encloser.part_id_),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)) {}
//=================================================================================================//
void TagBufferParticlesCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        buffer_particle_indicator_[index_i] = part_id_;
    }
}
//=================================================================================================//
template <class TargetPressure>
BufferInflowInjectionCK<TargetPressure>::
    BufferInflowInjectionCK(AlignedBoxPartByCell &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()), buffer_(buffer),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      sv_total_real_particles_(this->particles_->svTotalRealParticles()),
      create_real_particle_method_(this->particles_),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_buffer_particle_indicator_(this->particles_->template getVariableByName<int>("BufferParticleIndicator")),
      target_pressure_method_(this->particles_, this->sph_body_->getBaseMaterial()),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      upper_bound_fringe_(0.5 * this->sph_body_.getSPHBodyResolutionRef())
{
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
template <class TargetPressure>
template <class ExecutionPolicy, class EncloserType>
BufferInflowInjectionCK<TargetPressure>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : part_id_(encloser.part_id_),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      create_real_particle_(ex_policy, encloser.create_real_particle_method_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)),
      pressure_kernel_(ex_policy, encloser.target_pressure_method_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      upper_bound_fringe_(encloser.upper_bound_fringe_) {}
//=================================================================================================//
template <class TargetPressure>
void BufferInflowInjectionCK<TargetPressure>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (!aligned_box_->checkInBounds(pos_[index_i]))
    {
        if (aligned_box_->checkUpperBound(pos_[index_i], upper_bound_fringe_) &&
            buffer_particle_indicator_[index_i] == part_id_ &&
            index_i < *total_real_particles_)
        {
            UnsignedInt new_particle_index = create_real_particle_(index_i);
            buffer_particle_indicator_[new_particle_index] = 0;
            aligned_box_->imposeUpperPeriodic(pos_[index_i]);
            pressure_kernel_.reinitializeState(index_i, *physical_time_);
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BufferOutflowDeletionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      is_deltable_(encloser.part_id_, aligned_box_, pos_,
                   encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      remove_real_particle_(ex_policy, encloser.remove_real_particle_method_) {}
//=================================================================================================//
void BufferOutflowDeletionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (!aligned_box_->checkInBounds(pos_[index_i]))
    {
        if (is_deltable_(index_i) && index_i < *total_real_particles_)
        {
            remove_real_particle_(index_i, is_deltable_);
        }
    }
}
//=================================================================================================//
template <class TargetVelocity>
VelocityConditionCK<TargetVelocity>::
    VelocityConditionCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      target_velocity_method_(this->particles_),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")) {}
//=================================================================================================//
template <class TargetVelocity>
template <class ExecutionPolicy, class EncloserType>
VelocityConditionCK<TargetVelocity>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      transform_(&aligned_box_->getTransform()),
      velocity_kernel_(encloser.target_velocity_method_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class TargetVelocity>
void VelocityConditionCK<TargetVelocity>::UpdateKernel::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        Vecd frame_position = transform_->shiftBaseStationToFrame(pos_[index_i]);
        Vecd frame_velocity = transform_->xformBaseVecToFrame(vel_[index_i]);
        Vecd relaxed_frame_velocity = velocity_kernel_(frame_position, frame_velocity, *physical_time_);
        vel_[index_i] = transform_->xformFrameVecToBase(relaxed_frame_velocity);
    }
}
//=================================================================================================//
template <class KernelCorrectionType, class TargetPressure>
PressureConditionCK<KernelCorrectionType, TargetPressure>::
    PressureConditionCK(AlignedBoxPartByCell &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      ForcePriorCK(this->particles_, "BoundaryPressureForce"),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      kernel_correction_(this->particles_),
      target_pressure_method_(this->particles_),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_zero_gradient_residue_(this->particles_->template getVariableByName<Vecd>("ZeroGradientResidue")),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")) {}
//=================================================================================================//
template <class KernelCorrectionType, class TargetPressure>
template <class ExecutionPolicy, class EncloserType>
PressureConditionCK<KernelCorrectionType, TargetPressure>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : ForcePriorCK::UpdateKernel(ex_policy, encloser),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      correction_(ex_policy, encloser.kernel_correction_),
      pressure_kernel_(ex_policy, encloser.target_pressure_method_),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, class TargetPressure>
PressureConditionCK<KernelCorrectionType, TargetPressure>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    Vecd current_force = Vecd::Zero();
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        Vecd corrected_gradient_residue = correction_(index_i) * zero_gradient_residue_[index_i];
        current_force = pressure_kernel_.getTargetPressure(index_i, *physical_time_) *
                        mass_[index_i] * corrected_gradient_residue;
    }
    this->current_force_[index_i] = current_force;
    ForcePriorCK::UpdateKernel::update(index_i, dt);
}
//=================================================================================================//
template <typename ExecutionPolicy, class KernelCorrectionType, class TargetPressure, class TargetVelocity>
BidirectionalBoundaryCK<ExecutionPolicy, KernelCorrectionType, TargetPressure, TargetVelocity>::
    BidirectionalBoundaryCK(AlignedBoxPartByCell &emitter_by_cell, ParticleBuffer<Base> &inlet_buffer)
    : tag_buffer_particles_(emitter_by_cell),
      pressure_condition_(emitter_by_cell),
      emitter_injection_(emitter_by_cell, inlet_buffer),
      disposer_outflow_deletion_(emitter_by_cell) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // PRESSURE_BOUNDARY_CK_HPP
