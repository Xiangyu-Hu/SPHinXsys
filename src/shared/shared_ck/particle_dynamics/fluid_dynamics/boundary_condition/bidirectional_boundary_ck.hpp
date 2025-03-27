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
template <class BoundaryConditionType>
template <typename... Args>
BufferInflowInjectionCK<BoundaryConditionType>::
    BufferInflowInjectionCK(AlignedBoxPartByCell &aligned_box_part,
                            ParticleBuffer<Base> &buffer, Args &&...args)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()), buffer_(buffer),
      fluid_(DynamicCast<FluidType>(this, material)), condition_(std::forward<Args>(args)...),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      sv_total_real_particles_(this->particles_->svTotalRealParticles()),
      create_real_particle_method_(this->particles_),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_buffer_particle_indicator_(this->particles_->template getVariableByName<int>("BufferParticleIndicator")),
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
      recycle_condition_(ex_policy, encloser.boundary_condition_method_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      upper_bound_fringe_(encloser.upper_bound_fringe_) {}
//=================================================================================================//
template <class TargetPressure>
void BufferInflowInjectionCK<TargetPressure>::UpdateKernel::update(size_t index_i, Real dt)
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
            p_[index_i] = condition_.getPressure(p_[index_i], *physical_time_);
            rho_[index_i] = eos_.DensityFromPressure(p_[index_i]);
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
template <class KernelCorrectionType, typename ConditionType>
template <typename... Args>
PressureVelocityCondition<KernelCorrectionType, ConditionType>::
    BoundaryConditionCK(AlignedBoxPartByCell &aligned_box_part, Args &&...args)
    : BaseLocalDynamics<AlignedBoxPartByCell>(aligned_box_part),
      BaseStateCondition(this->particles_),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      kernel_correction_method_(this->particles_),
      condition_(std::forward<Args>(args)...),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
template <class ExecutionPolicy, class EncloserType>
BoundaryConditionCK<KernelCorrectionTypeType, TargetPressure>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseStateCondition(ex_policy, encloser),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      correction_kernel_(ex_policy, encloser.kernel_correction_method_),
      boundary_condition_(ex_policy, encloser.boundary_condition_method_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
PressureVelocityCondition<BoundaryConditionType>::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        Vecd corrected_gradient_residue = kernel_correction_(index_i) * zero_gradient_residue_[index_i];
        vel_[index_i] += dt * getPressure(p_[index_i], time) / rho_[index_i] * corrected_gradient_residue;

        Vecd frame_velocity = Vecd::Zero();
        frame_velocity[axis_] = transform_->xformBaseVecToFrame(vel_[index_i])[axis_];
        Vecd frame_position = transform_->shiftBaseStationToFrame(pos_[index_i]);
        frame_velocity[axis_] = getAxisVelocity(frame_position, frame_velocity[axis_], time);
        vel_[index_i] = transform_->xformFrameVecToBase(frame_velocity);
    }
}
//=================================================================================================//
template <typename ExecutionPolicy, class BoundaryConditionType>
template <typename... Args>
BidirectionalBoundaryCK<ExecutionPolicy, BoundaryConditionType>::
    BidirectionalBoundaryCK(AlignedBoxPartByCell &aligned_box_part,
                            ParticleBuffer<Base> &particle_buffer, Args &&...args)
    : tag_buffer_particles_(aligned_box_part),
      boundary_condition_(aligned_box_part, std::forward<Args>(args)...),
      inflow_injection_(aligned_box_part, particle_buffer, std::forward<Args>(args)...),
      outflow_deletion_(aligned_box_part) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // PRESSURE_BOUNDARY_CK_HPP
