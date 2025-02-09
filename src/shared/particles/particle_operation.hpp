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
SpawnRealParticle::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : total_real_particles_(encloser.sv_total_real_particles_->DelegatedData(ex_policy)),
      real_particles_bound_(encloser.real_particles_bound_),
      original_id_(encloser.dv_original_id_->DelegatedData(ex_policy))
{
    static_assert(std::is_base_of<SequencedPolicy, ExecutionPolicy>::value,
                  "SequencedPolicy is not the base of ExecutionPolicy!");
    OperationBetweenDataAssembles<ParticleVariables, DiscreteVariableArrays, DiscreteVariableArraysInitialization>
        initialize_discrete_variable_array;
    initialize_discrete_variable_array(encloser.variables_to_sort_, encloser.copyable_states_);
    OperationBetweenDataAssembles<DiscreteVariableArrays, VariableDataArrays, VariableDataArraysInitialization>
        initialize_variable_data_array;
    initialize_variable_data_array(encloser.copyable_states_, copyable_state_data_arrays_, ex_policy);
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_OPERATION_HPP