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
      ce_mesh_(createConstantEntity<Mesh>(execution_policy, "Mesh", DynamicCast<CellLinkedListType>(this, real_body.getCellLinkedList()))),
      sv_total_real_particles_(particles_->getSingularVariableByName<UnsignedInt>("TotalRealParticles")),
      dv_particle_offset_list_(nullptr), mesh_(ce_mesh_->DataAddress(execution_policy)),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      particle_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleIDList", particles_->ParticlesBound())),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleOffsetList", mesh_->NumberOfCells() + 1)),
      current_size_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "CurrentCellSize", mesh_->NumberOfCells()))
{
    dv_particle_offset_list_ = particles_->getVariableByName<UnsignedInt>("ParticleOffsetList");
    particles_->addVariableToWrite<UnsignedInt>("ParticleIDList");
}
//=================================================================================================//
template <typename CellLinkedListType>
UnsignedInt *UpdateCellLinkedList<CellLinkedListType>::setParticleOffsetListUpperBound()
{
    UnsignedInt total_real_particles = particles_->TotalRealParticles();
    UnsignedInt *particle_offset_list = dv_particle_offset_list_->DataField();
    particle_offset_list[mesh_->NumberOfCells()] = total_real_particles;
    return particle_offset_list;
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
