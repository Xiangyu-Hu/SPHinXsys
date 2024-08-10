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
    : LocalDynamics(real_body),
      mesh_(createConstantEntity<Mesh>(execution_policy, "Mesh", DynamicCast<CellLinkedListType>(this, real_body.getCellLinkedList()))),
      number_of_cells_(mesh_->NumberOfCells()), particles_bound_(particles_->ParticlesBound()),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      particle_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleIDList", particles_bound_)),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleOffsetList", number_of_cells_ + 1)),
      current_size_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "CurrentCellSize", number_of_cells_))
{
    particles_->addVariableToWrite<UnsignedInt>("ParticleIDList");
}
//=================================================================================================//
template <class CellLinkedListType, class ExecutionPolicy>
UpdateCellLinkedList<CellLinkedListType, ExecutionPolicy>::UpdateCellLinkedList(RealBody &real_body)
    : UpdateCellLinkedList<CellLinkedListType>(ExecutionPolicy{}, real_body),
      BaseDynamics<void>(){};
//=================================================================================================//
template <class CellLinkedListType, class ExecutionPolicy>
void UpdateCellLinkedList<CellLinkedListType, ExecutionPolicy>::exec(Real dt)
{
    this->clearParticleOffsetList(ExecutionPolicy{});
    this->incrementCellSize(ExecutionPolicy{});
    this->exclusiveScanParticleOffsetList(ExecutionPolicy{});
    this->updateCellLists(ExecutionPolicy{});
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
