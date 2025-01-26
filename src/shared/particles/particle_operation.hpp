#ifndef PARTICLE_OPERATION_HPP
#define PARTICLE_OPERATION_HPP

#include "particle_operation.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
void CopyParticleStateCK::operator()(
    VariableAllocationPair<AllocatedDataArray<DataType>> &variable_allocation_pair,
    size_t index, size_t another_index)
{
    for (size_t i = 0; i != variable_allocation_pair.second; ++i)
    {
        variable_allocation_pair.first[i][index] = variable_allocation_pair.first[i][another_index];
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
CreateRealParticleFrom::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      real_particles_bound_(encloser.real_particles_bound_),
      original_id_(encloser.dv_original_id_->DelegatedData(ex_policy))//,
      copy_particle_state_(
        this->initializeCopyParticleState(ex_policy, encloser.variables_to_sort_, encloser.copyable_states_)) {}
//=================================================================================================//
template <class ExecutionPolicy>
OperationOnDataAssemble<VariableDataArrays, CopyParticleStateCK>
CreateRealParticleFrom::ComputingKernel::initializeCopyParticleState(
    const ExecutionPolicy &ex_policy,
    ParticleVariables &variables_to_sort, DiscreteVariableArrays &copyable_states)
{
    OperationBetweenDataAssembles<ParticleVariables, DiscreteVariableArrays, DiscreteVariableArraysInitialization>
        initialize_discrete_variable_array(variables_to_sort, copyable_states);
    initialize_discrete_variable_array();
    OperationBetweenDataAssembles<DiscreteVariableArrays, VariableDataArrays, VariableDataArraysInitialization>
        initialize_variable_data_array(copyable_states, copyable_state_data_arrays_);
    initialize_variable_data_array(ex_policy);
    return OperationOnDataAssemble<VariableDataArrays, CopyParticleStateCK>(copyable_state_data_arrays_);
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_OPERATION_HPP