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
    : LocalDynamics(real_body), BaseDynamics<void>(real_body),
      mesh_(real_body.getCellLinkedList()),
      number_of_cells_(mesh_.NumberOfCells()),
      particles_bound_(particles_->ParticlesBound()),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      particle_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, particles_bound_, "ParticleIDList")),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, number_of_cells_ + 1, "ParticleOffsetList")),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, number_of_cells_, "CurrentCellSize")),
      computing_kernel_buffer_(nullptr) {}
//=================================================================================================//
template <typename CellLinkedListType>
UpdateCellLinkedList<CellLinkedListType>::ComputingKernel::
    ComputingKernel(UpdateCellLinkedList<CellLinkedListType> &update_cell_linked_list)
    : mesh_(update_cell_linked_list.mesh_),
      number_of_cells_(update_cell_linked_list.number_of_cells_),
      particles_bound_(update_cell_linked_list.particles_bound_),
      pos_(update_cell_linked_list.pos_),
      particle_id_list_(update_cell_linked_list.particle_id_list_),
      particle_offset_list_(update_cell_linked_list.particle_offset_list_) {}
//=================================================================================================//
template <class CellLinkedListType, class ExecutionPolicy>
UpdateCellLinkedList<CellLinkedListType, ExecutionPolicy>::UpdateCellLinkedList(RealBody &real_body)
    : UpdateCellLinkedList<CellLinkedListType>(ExecutionPolicy{}, real_body),
      BaseDynamics<void>(), kernel_implementation_(*this){};
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
