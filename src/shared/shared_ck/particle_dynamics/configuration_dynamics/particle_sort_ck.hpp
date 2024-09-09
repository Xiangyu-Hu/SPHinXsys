#ifndef PARTICLE_SORT_HPP
#define PARTICLE_SORT_HPP

#include "particle_sort_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
void UpdateSortableVariables::InitializeTemporaryVariables::operator()(
    UniquePtrKeeper<DiscreteVariable<DataType>> &variable_ptr_keeper, UnsignedInt data_size)
{
    variable_ptr_keeper.template createPtr<DiscreteVariable<DataType>>("Temporary", data_size);
}

//=================================================================================================//
template <class ExecutionPolicy, typename DataType>
void UpdateSortableVariables::operator()(
    DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
    ExecutionPolicy &ex_policy, DiscreteVariable<UnsignedInt> *dv_index_permutation)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    DataType *temp_data_field = std::get<type_index>(temp_variables_)->DelegatedDataField(ex_policy);

    UnsignedInt *index_permutation = dv_index_permutation->DelegatedDataField(ex_policy);

    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    for (size_t k = 0; k != variables.size(); ++k)
    {
        DataType *sorted_data_field = variables[k]->DelegatedDataField(ex_policy);
        particle_for(ex_policy, IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { temp_data_field[i] = sorted_data_field[i]; });
        particle_for(ex_policy, IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { sorted_data_field[i] = temp_data_field[index_permutation[i]]; });
    }
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SORT_HPP