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
EmitterInflowInjectionCK::
    EmitterInflowInjectionCK(AlignedBoxPartByParticle &aligned_box_part, ParticleBuffer<Base> &buffer)
    : BaseLocalDynamics<AlignedBoxPartByParticle>(aligned_box_part),
      buffer_(buffer), sv_aligned_box_(aligned_box_part.svAlignedBox()),
      rho0_(particles_->getBaseMaterial().ReferenceDensity()),
      copyable_states_(),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_rho_(particles_->getVariableByName<Real>("Density")),
      dv_p_(particles_->getVariableByName<Real>("Pressure"))
{
    OperationBetweenDataAssembles<ParticleVariables, DiscreteVariableArrays, DiscreteVariableArraysInitialization>
        initialize_discrete_variable_array(particles_->VariablesToSort(), copyable_states_);
    initialize_discrete_variable_array();
    buffer_.checkParticlesReserved();
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
EmitterInflowInjectionCK::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
      rho0_(encloser.rho_), copyable_state_data_arrays_(),
      copy_particle_state_(copyable_state_data_arrays_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      p_(encloser.dv_p_->DelegatedData(ex_policy)),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy))
{
    OperationBetweenDataAssembles<DiscreteVariableArrays, VariableDataArrays, VariableDataArraysInitialization>
        initialize_variable_data_array(encloser.copyable_states_, copyable_state_data_arrays_);
    initialize_variable_data_array();
}
//=================================================================================================//
void EmitterInflowInjectionCK::UpdateKernel::update(size_t index_i, Real dt)
{
    if (aligned_box_->checkUpperBound(pos_[index_i]))
    {
        particles_->createRealParticleFrom(index_i);

        /** Periodic bounding. */
        pos_[index_i] = aligned_box_.getUpperPeriodic(pos_[index_i]);
        rho_[index_i] = fluid_.ReferenceDensity();
        p_[index_i] = fluid_.getPressure(rho_[index_i]);
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_HPP
