#ifndef BIDIRECTIONAL_BOUNDARY_CK_HPP
#define BIDIRECTIONAL_BOUNDARY_CK_HPP

#include "bidirectional_boundary_ck.h"

#include "complex_geometry.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BufferIndicationCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : part_id_(encloser.part_id_),
      oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_indicator_(encloser.dv_buffer_indicator_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline void BufferIndicationCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (oriented_box_->checkContain(pos_[index_i]))
    {
        buffer_indicator_[index_i] = part_id_;
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
ResetBufferCorrectionMatrixCK::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      B_(encloser.dv_B_->DelegatedData(ex_policy)),
      radius_(encloser.radius_) {}
//=================================================================================================//
inline void ResetBufferCorrectionMatrixCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (oriented_box_->checkLowerBound(pos_[index_i], -radius_))
    {
        B_[index_i] = Matd::Identity();
    }
}
//=================================================================================================//
template <class ConditionType>
template <typename... Args>
BufferInflowInjectionCK<ConditionType>::
    BufferInflowInjectionCK(OrientedBoxByCell &oriented_box_part, Args &&...args)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      part_id_(oriented_box_part.getPartID()),
      fluid_(DynamicCast<FluidType>(this, sph_body_->getMatterMaterial())),
      condition_(std::forward<Args>(args)...),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
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
    : part_id_(encloser.part_id_), eos_(ex_policy, encloser.fluid_), condition_(encloser.condition_),
      oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
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
    if (!oriented_box_->checkInBounds(pos_[index_i]))
    {
        if (oriented_box_->checkUpperBound(pos_[index_i], upper_bound_fringe_) &&
            buffer_indicator_[index_i] == part_id_ &&
            index_i < *total_real_particles_)
        {
            UnsignedInt new_particle_index = spawn_real_particle_(index_i);
            buffer_indicator_[new_particle_index] = 0;
            pos_[index_i] = oriented_box_->getUpperPeriodic(pos_[index_i]);
            p_[index_i] = condition_.getPressure(p_[index_i], *physical_time_);
            rho_[index_i] = eos_.DensityFromPressure(index_i, p_[index_i]);
        }
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
BufferOutflowIndication::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      life_status_(encloser.dv_life_status_->DelegatedData(ex_policy)),
      is_deltable_(encloser.part_id_, oriented_box_, pos_,
                   encloser.dv_buffer_indicator_->DelegatedData(ex_policy)),
      total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline bool BufferOutflowIndication::UpdateKernel::IsDeletable::operator()(size_t index_i) const
{
    return buffer_indicator_[index_i] == part_id_ &&
           oriented_box_->checkLowerBound(pos_[index_i]);
}
//=================================================================================================//
inline void BufferOutflowIndication::UpdateKernel::update(size_t index_i, Real dt)
{
    if (!oriented_box_->checkInBounds(pos_[index_i]))
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
    PressureVelocityCondition(OrientedBoxByCell &oriented_box_part, Args &&...args)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      BaseStateCondition(this->particles_),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
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
      oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      correction_kernel_(ex_policy, encloser.kernel_correction_method_),
      condition_(encloser.condition_),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      kernel_gradient_integral_(encloser.dv_kernel_gradient_integral_->DelegatedData(ex_policy)),
      axis_(oriented_box_->ReferenceAxis()), transform_(&oriented_box_->getTransform()) {}
//=================================================================================================//
template <class KernelCorrectionType, typename ConditionType>
void PressureVelocityCondition<KernelCorrectionType, ConditionType>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (oriented_box_->checkContain(pos_[index_i]))
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
template <typename ConditionType>
template <typename... Args>
SupplementaryCondition<ConditionType>::SupplementaryCondition(
    OrientedBoxByCell &oriented_box_part, Args &&...args)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      condition_method_(this->particles_, std::forward<Args>(args)...),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")) {}
//=================================================================================================//
template <typename ConditionType>
template <class ExecutionPolicy, class EncloserType>
SupplementaryCondition<ConditionType>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      condition_(ex_policy, encloser.condition_method_),
      pos_(encloser.dv_pos_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename ConditionType>
void SupplementaryCondition<ConditionType>::UpdateKernel::update(size_t index_i, Real dt)
{
    if (oriented_box_->checkContain(pos_[index_i]))
    {
        condition_(index_i, dt);
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class ConditionType, typename... Args>
AbstractBidirectionalBoundary &AbstractBidirectionalBoundary::
    addSupplementaryCondition(OrientedBoxByCell &oriented_box_part, Args &&...args)
{
    auto *condition = supplementary_conditions_keeper_.template createPtr<
        StateDynamics<ExecutionPolicy, SupplementaryCondition<ConditionType>>>(
        oriented_box_part, std::forward<Args>(args)...);
    supplementary_conditions_.push_back(condition);
    return *this;
}
//=================================================================================================//
template <typename ExecutionPolicy, class KernelCorrectionType, class ConditionType>
template <typename... Args>
BidirectionalBoundaryCK<ExecutionPolicy, KernelCorrectionType, ConditionType>::
    BidirectionalBoundaryCK(OrientedBoxByCell &oriented_box_part, Args &&...args)
    : AbstractBidirectionalBoundary(), tag_buffer_particles_(oriented_box_part),
      reset_buffer_correction_matrix_(oriented_box_part),
      boundary_condition_(oriented_box_part, std::forward<Args>(args)...),
      inflow_injection_(oriented_box_part, std::forward<Args>(args)...),
      outflow_indication_(oriented_box_part) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // BIDIRECTIONAL_BOUNDARY_CK_HPP
