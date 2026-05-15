#ifndef PARTICLE_OPERATION_HPP
#define PARTICLE_OPERATION_HPP

#include "particle_operation.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
void CopyParticleStateCK::operator()(
    DataArray<DataType> &variable_data_array, size_t index, size_t another_index)
{
    for (size_t i = 0; i != variable_data_array.ArraySize(); ++i)
    {
        variable_data_array[i][index] = variable_data_array[i][another_index];
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
SpawnRealParticle::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      particles_bound_(encloser.particles_bound_),
      original_id_(encloser.dv_original_id_->DelegatedData(ex_policy))
{
    OperationBetweenDataAssembles<DiscreteVariables, VariableArrayAssemble, VariableArrayAssembleInitialization>
        initialize_discrete_variable_array;
    initialize_discrete_variable_array(encloser.evolving_variables_, encloser.copyable_states_);
    OperationBetweenDataAssembles<VariableArrayAssemble, VariableDataArrayAssemble, VariableDataArrayAssembleInitialization>
        initialize_variable_data_array;
    initialize_variable_data_array(encloser.copyable_states_, copyable_state_data_arrays_, ex_policy);
}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
RemoveRealParticle::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      original_id_(encloser.dv_original_id_->DelegatedData(ex_policy))
{
    OperationBetweenDataAssembles<DiscreteVariables, VariableArrayAssemble, VariableArrayAssembleInitialization>
        initialize_discrete_variable_array;
    initialize_discrete_variable_array(encloser.evolving_variables_, encloser.copyable_states_);
    OperationBetweenDataAssembles<VariableArrayAssemble, VariableDataArrayAssemble, VariableDataArrayAssembleInitialization>
        initialize_variable_data_array;
    initialize_variable_data_array(encloser.copyable_states_, copyable_state_data_arrays_, ex_policy);
}
} // namespace SPH
#endif // PARTICLE_OPERATION_HPP