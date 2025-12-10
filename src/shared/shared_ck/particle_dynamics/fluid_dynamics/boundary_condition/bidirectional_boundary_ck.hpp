#ifndef BIDIRECTIONAL_BOUNDARY_CK_HPP
#define BIDIRECTIONAL_BOUNDARY_CK_HPP

#include "bidirectional_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BufferIndicationCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : part_id_(encloser.part_id_),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_indicator_(encloser.dv_buffer_indicator_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void BufferIndicationCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        buffer_indicator_[index_i] = part_id_;
    }
}
//=================================================================================================//
template <class ConditionType>
template <typename... Args>
BufferInflowInjectionCK<ConditionType>::
    BufferInflowInjectionCK(AlignedBoxByCell &aligned_box_part, Args &&...args)
    : BaseLocalDynamics<AlignedBoxByCell>(aligned_box_part),
      part_id_(aligned_box_part.getPartID()),
      fluid_(DynamicCast<FluidType>(this, sph_body_->getBaseMaterial())),
      condition_(std::forward<Args>(args)...),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      sv_total_real_particles_(this->particles_->svTotalRealParticles()),
      spawn_real_particle_method_(this->particles_),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      dv_buffer_indicator_(this->particles_->template getVariableByName<int>("BufferIndicator")),
      sv_physical_time_(this->sph_system_->template getSystemVariableByName<Real>("PhysicalTime")),
      dv_p_(this->particles_->template getVariableByName<Real>("Pressure")),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      upper_bound_fringe_(0.5 * this->sph_body_->getSPHBodyResolutionRef())
{
    particles_->checkEnoughReserve();
}
//=================================================================================================//
template <class ConditionType>
template <class ExecutionPolicy, class EncloserType>
BufferInflowInjectionCK<ConditionType>::
    UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : part_id_(encloser.part_id_), eos_(encloser.fluid_), condition_(encloser.condition_),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      spawn_real_particle_(ex_policy, encloser.spawn_real_particle_method_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_indicator_(encloser.dv_buffer_indicator_->DelegatedData(ex_policy)),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      upper_bound_fringe_(encloser.upper_bound_fringe_) {}
//=================================================================================================//
template <class ConditionType>
void BufferInflowInjectionCK<ConditionType>::UpdateKernel::update(size_t index_i, Real dt)
{
    if (!aligned_box_->checkInBounds(pos_[index_i]))
    {
        if (aligned_box_->checkUpperBound(pos_[index_i], upper_bound_fringe_) &&
            buffer_indicator_[index_i] == part_id_ &&
            index_i < *total_real_particles_)
        {
            UnsignedInt new_particle_index = spawn_real_particle_(index_i);
            buffer_indicator_[new_particle_index] = 0;
            pos_[index_i] = aligned_box_->getUpperPeriodic(pos_[index_i]);
            p_[index_i] = condition_.getPressure(p_[index_i], *physical_time_);
            rho_[index_i] = eos_.DensityFromPressure(p_[index_i]);
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BufferOutflowIndication::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      life_status_(encloser.dv_life_status_->DelegatedData(ex_policy)),
      is_deltable_(encloser.part_id_, aligned_box_, pos_,
                   encloser.dv_buffer_indicator_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void BufferOutflowIndication::UpdateKernel::update(size_t index_i, Real dt)
{
    if (!aligned_box_->checkInBounds(pos_[index_i]))
    {
        if (is_deltable_(index_i) && index_i < *total_real_particles_)
        {
            life_status_[index_i] = 1; // mark as to delete but will not delete immediately
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
OutflowParticleDeletion::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : remove_real_particle_(ex_policy, encloser.remove_real_particle_method_),
      life_status_(encloser.dv_life_status_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void OutflowParticleDeletion::UpdateKernel::update(UnsignedInt index_i, Real dt)
{
    if (life_status_[index_i] == 1) // to delete
    {
        remove_real_particle_(index_i, life_status_);
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
template <typename... Args>
PressureVelocityCondition<KernelCorrectionType, ConditionType>::
    PressureVelocityCondition(AlignedBoxByCell &aligned_box_part, Args &&...args)
    : BaseLocalDynamics<AlignedBoxByCell>(aligned_box_part),
      BaseStateCondition(this->particles_),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      kernel_correction_method_(this->particles_),
      condition_(std::forward<Args>(args)...),
      sv_physical_time_(this->sph_system_->template getSystemVariableByName<Real>("PhysicalTime")),
      dv_kernel_gradient_integral_(this->particles_->template getVariableByName<Vecd>("KernelGradientIntegral")) {}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
template <class ExecutionPolicy, class EncloserType>
PressureVelocityCondition<KernelCorrectionType, ConditionType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseStateCondition::ComputingKernel(ex_policy, encloser),
      aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      correction_kernel_(ex_policy, encloser.kernel_correction_method_),
      condition_(encloser.condition_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      kernel_gradient_integral_(encloser.dv_kernel_gradient_integral_->DelegatedData(ex_policy)),
      axis_(aligned_box_->AlignmentAxis()), transform_(&aligned_box_->getTransform()) {}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
void PressureVelocityCondition<KernelCorrectionType, ConditionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        Vecd corrected_residual = correction_kernel_(index_i) * kernel_gradient_integral_[index_i];
        vel_[index_i] -= dt * condition_.getPressure(p_[index_i], *physical_time_) /
                         rho_[index_i] * corrected_residual;

        Vecd frame_velocity = Vecd::Zero();
        frame_velocity[axis_] = transform_->xformBaseVecToFrame(vel_[index_i])[axis_];
        Vecd frame_position = transform_->shiftBaseStationToFrame(pos_[index_i]);
        frame_velocity[axis_] = condition_.getAxisVelocity(frame_position, frame_velocity[axis_], *physical_time_);
        vel_[index_i] = transform_->xformFrameVecToBase(frame_velocity);
    }
}
//=================================================================================================//
template <typename ExecutionPolicy, class KernelCorrectionType, class ConditionType>
template <typename... Args>
BidirectionalBoundaryCK<ExecutionPolicy, KernelCorrectionType, ConditionType>::
    BidirectionalBoundaryCK(AlignedBoxByCell &aligned_box_part, Args &&...args)
    : AbstractDynamics(), tag_buffer_particles_(aligned_box_part),
      boundary_condition_(aligned_box_part, std::forward<Args>(args)...),
      inflow_injection_(aligned_box_part, std::forward<Args>(args)...),
      outflow_indication_(aligned_box_part) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BOUNDARY_CK_HPP
