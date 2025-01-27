#ifndef PARTICLE_SORT_HPP
#define PARTICLE_SORT_HPP

#include "particle_sort_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
void UpdateSortableVariables::InitializeTemporaryVariables::operator()(
    UniquePtr<DiscreteVariable<DataType>> &variable_ptr, UnsignedInt data_size)
{
    variable_ptr = makeUnique<DiscreteVariable<DataType>>("Temporary", data_size);
}
//=================================================================================================//
template <class ExecutionPolicy, typename DataType>
void UpdateSortableVariables::operator()(
    DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables,
    ExecutionPolicy &ex_policy, UnsignedInt total_real_particles,
    DiscreteVariable<UnsignedInt> *dv_index_permutation)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    DataType *temp_data_field = std::get<type_index>(temp_variables_)->DelegatedData(ex_policy);

    UnsignedInt *index_permutation = dv_index_permutation->DelegatedData(ex_policy);

    for (size_t k = 0; k != variables.size(); ++k)
    {
        DataType *sorted_data_field = variables[k]->DelegatedData(ex_policy);
        particle_for(ex_policy, IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { temp_data_field[i] = sorted_data_field[i]; });
        particle_for(ex_policy, IndexRange(0, total_real_particles),
                     [=](size_t i)
                     { sorted_data_field[i] = temp_data_field[index_permutation[i]]; });
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
QuickSort::QuickSort(const ExecutionPolicy &ex_policy,
                     DiscreteVariable<UnsignedInt> *dv_sequence,
                     DiscreteVariable<UnsignedInt> *dv_index_permutation)
    : sequence_(dv_sequence->DelegatedData(ex_policy)),
      index_permutation_(dv_index_permutation->DelegatedData(ex_policy)),
      swap_particle_index_(sequence_, index_permutation_), compare_(),
      quick_sort_particle_range_(sequence_, 0, compare_, swap_particle_index_),
      quick_sort_particle_body_() {}
//=================================================================================================//
template <class ExecutionPolicy, class SortMethodType>
ParticleSortCK<ExecutionPolicy, SortMethodType>::ParticleSortCK(RealBody &real_body)
    : LocalDynamics(real_body), BaseDynamics<void>(),
      ex_policy_(ExecutionPolicy{}),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())),
      mesh_(cell_linked_list_),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_sequence_(particles_->registerDiscreteVariableOnly<UnsignedInt>(
          "Sequence", particles_->ParticlesBound())),
      dv_index_permutation_(particles_->registerDiscreteVariableOnly<UnsignedInt>(
          "IndexPermutation", particles_->ParticlesBound())),
      dv_original_id_(particles_->getVariableByName<UnsignedInt>("OriginalID")),
      dv_sorted_id_(particles_->getVariableByName<UnsignedInt>("SortedID")),
      update_variables_to_sort_(particles_),
      sort_method_(ExecutionPolicy{}, dv_sequence_, dv_index_permutation_),
      kernel_implementation_(*this)
{
    particles_->addVariableToSort<UnsignedInt>("OriginalID");
}
//=================================================================================================//
template <class ExecutionPolicy, class SortMethodType>
ParticleSortCK<ExecutionPolicy, SortMethodType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    ParticleSortCK<ExecutionPolicy, SortMethodType> &encloser)
    : mesh_(encloser.mesh_), pos_(encloser.dv_pos_->DelegatedData(ex_policy)),
      sequence_(encloser.dv_sequence_->DelegatedData(ex_policy)),
      index_permutation_(encloser.dv_index_permutation_->DelegatedData(ex_policy)),
      original_id_(encloser.dv_original_id_->DelegatedData(ex_policy)),
      sorted_id_(encloser.dv_sorted_id_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class SortMethodType>
void ParticleSortCK<ExecutionPolicy, SortMethodType>::ComputingKernel::
    prepareSequence(UnsignedInt index_i)
{
    sequence_[index_i] = mesh_.transferMeshIndexToMortonOrder(mesh_.CellIndexFromPosition(pos_[index_i]));
    index_permutation_[index_i] = index_i;
}
//=================================================================================================//
template <class ExecutionPolicy, class SortMethodType>
void ParticleSortCK<ExecutionPolicy, SortMethodType>::ComputingKernel::
    updateSortedID(UnsignedInt index_i)
{
    sorted_id_[original_id_[index_i]] = index_i;
}
//=================================================================================================//
template <class ExecutionPolicy, class SortMethodType>
void ParticleSortCK<ExecutionPolicy, SortMethodType>::exec(Real dt)
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();

    particle_for(ex_policy_, IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->prepareSequence(i); });

    sort_method_.sort(ex_policy_, particles_);
    update_variables_to_sort_(particles_->VariablesToSort(), ex_policy_, total_real_particles, dv_index_permutation_);

    particle_for(ex_policy_, IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateSortedID(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SORT_HPP