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
template <class ExecutionPolicy, class EncloserType>
DisposerOutflowDeletionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      remove_real_particle_(ex_policy, encloser.remove_real_particle_method_),
      rho0_(encloser.rho0_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)) {}
//=================================================================================================//
void DisposerOutflowDeletionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        remove_real_particle_(index_i);
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TagBufferParticlesCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      //   part_id_(encloser.sv_part_id_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy))
{
}
//=================================================================================================//
void TagBufferParticlesCK::UpdateKernel::update(size_t index_i, Real dt)
{
    int buffer_indicator = 0;
    if (aligned_box_->checkContain(pos_[index_i]))
    {
        buffer_indicator = 1;
    }
    // TODO: we need a getPartID sv version!
    // TODO: maybe a initial state to set all buffer indicator = 0 first?
    buffer_particle_indicator_[index_i] = buffer_indicator;
}
//=================================================================================================//
template <class AlignedBoxPartType, class KernelCorrectionType, class ConditionFunction>
PressureConditionCK<AlignedBoxPartType, KernelCorrectionType, ConditionFunction>::
    PressureConditionCK(AlignedBoxPartType &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      sv_aligned_box_(aligned_box_part.svAlignedBox()),
      condition_function_(this->particles_),
      dv_buffer_particle_indicator_(this->particles_->template getVariableByName<int>("BufferParticleIndicator")),
      dv_zero_gradient_residue_(this->particles_->template getVariableByName<Vecd>("ZeroGradientResidue")),
      dv_pos_(this->particles_->template getVariableByName<Vecd>("Position")),
      sv_physical_time_(this->sph_system_.template getSystemVariableByName<Real>("PhysicalTime")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      kernel_correction_(this->particles_),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density"))
{
}

//=================================================================================================//
template <class AlignedBoxPartType, class KernelCorrectionType, class ConditionFunction>
template <class ExecutionPolicy, class EncloserType>
PressureConditionCK<AlignedBoxPartType, KernelCorrectionType, ConditionFunction>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      condition_(ex_policy, encloser.condition_function_),
      buffer_particle_indicator_(encloser.dv_buffer_particle_indicator_->DelegatedData(ex_policy)),
      zero_gradient_residue_(encloser.dv_zero_gradient_residue_->DelegatedData(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      physical_time_(encloser.sv_physical_time_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)),
      correction_(ex_policy, encloser.kernel_correction_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      transform_(&aligned_box_->getTransform())
{
}
//=================================================================================================//
template <class AlignedBoxPartType, class KernelCorrectionType, class ConditionFunction>
void PressureConditionCK<AlignedBoxPartType, KernelCorrectionType, ConditionFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (buffer_particle_indicator_[index_i] != 0)
    {
        if (aligned_box_->checkContain(pos_[index_i]))
        {
            // vel_[index_i] += 2.0 * kernel_sum_[index_i] * target_pressure_(p_[index_i], *physical_time_) / rho_[index_i] * dt;
            Real pressure = condition_(index_i, *physical_time_);
            vel_[index_i] += correction_(index_i) * this->zero_gradient_residue_[index_i] * pressure / rho_[index_i] * dt;

            Vecd frame_velocity = Vecd::Zero();

            frame_velocity[xAxis] = this->transform_->xformBaseVecToFrame(vel_[index_i])[xAxis];
            vel_[index_i] = transform_->xformFrameVecToBase(frame_velocity);
        }
    }
}
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_HPP
