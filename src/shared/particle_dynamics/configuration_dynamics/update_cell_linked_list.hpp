#ifndef UPDATE_CELL_LINKED_LIST_HPP
#define UPDATE_CELL_LINKED_LIST_HPP

#include "update_cell_linked_list.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
template <typename CellLinkedListType>
template <class ExecutionPolicy>
UpdateCellLinkedList<CellLinkedListType>::
    UpdateCellLinkedList(const ExecutionPolicy &execution_policy, RealBody &real_body)
    : LocalDynamics(real_body), mesh_(DynamicCast<CellLinkedListType>(this, real_body.getCellLinkedList())),
      number_of_cells_(mesh_.NumberOfCells()), particles_bound_(particles_->ParticlesBound()),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      particle_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleIDList", particles_bound_)),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleOffsetList", number_of_cells_ + 1)),
      current_size_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "CurrentCellSize", number_of_cells_)),
      v_total_real_particles_(particles_->getSingularVariableByName<UnsignedInt>("TotalRealParticles")),
      v_particle_offset_list_(particles_->getVariableByName<UnsignedInt>("ParticleOffsetList")) {}
//=================================================================================================//
template <typename CellLinkedListType>
UnsignedInt *UpdateCellLinkedList<CellLinkedListType>::setParticleOffsetListUpperBound()
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    UnsignedInt *particle_offset_list = v_particle_offset_list_->DataField();
    particle_offset_list[number_of_cells_] = total_real_particles;
    return particle_offset_list;
}
//=================================================================================================//
template <typename CellLinkedListType>
template <class ExecutionPolicy>
void UpdateCellLinkedList<CellLinkedListType>::
    setParticleOffsetListUpperBound(const ExecutionPolicy &execution_policy)
{
    setParticleOffsetListUpperBound();
}
//=================================================================================================//
template <typename CellLinkedListType>
UpdateCellLinkedList<CellLinkedListType>::ComputingKernel::
    ComputingKernel(UpdateCellLinkedList<CellLinkedListType> &update_cell_linked_list)
    : mesh_(update_cell_linked_list.mesh_),
      number_of_cells_(update_cell_linked_list.number_of_cells_),
      particles_bound_(update_cell_linked_list.particles_bound_),
      pos_(update_cell_linked_list.pos_),
      particle_id_list_(update_cell_linked_list.particle_id_list_),
      particle_offset_list_(update_cell_linked_list.particle_offset_list_),
      current_size_list_(update_cell_linked_list.current_size_list_),
      total_real_particles_(update_cell_linked_list.v_total_real_particles_->ValueAddress()) {}
//=================================================================================================//
template <typename CellLinkedListType>
void UpdateCellLinkedList<CellLinkedListType>::ComputingKernel::
    clearOffsetLists(UnsignedInt linear_cell_index)
{
    particle_offset_list_[linear_cell_index] = 0;
    current_size_list_[linear_cell_index] = 0;
}
//=================================================================================================//
template <typename CellLinkedListType>
void UpdateCellLinkedList<CellLinkedListType>::ComputingKernel::
    incrementCellSize(UnsignedInt particle_i)
{
    UnsignedInt linear_cell_index = mesh_.LinearCellIndexFromPosition(pos_[particle_i]);
    AtomicUnsignedIntRef<ExecutionPolicy>::type atomic_cell_size(offset_cell_size[linear_cell_index]);
    ++atomic_cell_size;
} //=================================================================================================//
template <typename CellLinkedListType>
void UpdateCellLinkedList<CellLinkedListType>::ComputingKernel::
    updateCellLists(UnsignedInt particle_i)
{
    UnsignedInt linear_cell_index = mesh_.LinearCellIndexFromPosition(pos_[particle_i]);
    AtomicUnsignedIntRef<ExecutionPolicy>::type atomic_current_cell_size(curr_cell_size[linear_cell_index]);
    particle_id_list[offset_cell_size[linear_cell_index] + atomic_current_cell_size++] = particle_i;
}
//=================================================================================================//
template <class CellLinkedListType, class ExecutionPolicy>
UpdateCellLinkedList<CellLinkedListType, ExecutionPolicy>::UpdateCellLinkedList(RealBody &real_body)
    : UpdateCellLinkedList<CellLinkedListType>(ExecutionPolicy{}, real_body),
      BaseDynamics<void>(), kernel_implementation_(*this){};
//=================================================================================================//
template <class CellLinkedListType, class ExecutionPolicy>
void UpdateCellLinkedList<CellLinkedListType, ExecutionPolicy>::exec(Real dt)
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    ComputingKernel *computing_kernel = kernel_implementation_.getComputingKernel();

    setParticleOffsetListUpperBound(ExecutionPolicy{});
    particle_for(ExecutionPolicy{},
                 IndexRange(0, number_of_cells_),
                 [=](size_t i)
                 { computing_kernel->clearOffsetLists(i); });

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->incrementCellSize(i); });

    computing_kernel->exclusiveScanParticleOffsetList(ExecutionPolicy{});

    particle_for(ExecutionPolicy{},
                 IndexRange(0, total_real_particles),
                 [=](size_t i)
                 { computing_kernel->updateCellLists(i); });
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
