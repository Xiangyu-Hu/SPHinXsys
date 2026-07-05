#ifndef PARTICLE_OPERATION_HPP
#define PARTICLE_OPERATION_HPP

#include "particle_operation.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
void CopyParticleStateCK::operator()(
    VariableArrayView<DataType> &variable_array_view, UnsignedInt index, UnsignedInt another_index)
{
    for (UnsignedInt i = 0; i != variable_array_view.ArraySize(); ++i)
    {
        for (UnsignedInt j = 0; j != variable_array_view[i].Width(); ++j)
        {
            variable_array_view[i][index][j] = variable_array_view[i][another_index][j];
        }
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
    OperationBetweenDataAssembles<VariableArrayAssemble, VariableArrayViewAssemble, VariableArrayViewAssembleInitialization>
        initialize_variable_array_view;
    initialize_variable_array_view(encloser.copyable_states_, copyable_state_data_arrays_, ex_policy);
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
    OperationBetweenDataAssembles<VariableArrayAssemble, VariableArrayViewAssemble, VariableArrayViewAssembleInitialization>
        initialize_variable_array_view;
    initialize_variable_array_view(encloser.copyable_states_, copyable_state_data_arrays_, ex_policy);
}
} // namespace SPH
#endif // PARTICLE_OPERATION_HPP