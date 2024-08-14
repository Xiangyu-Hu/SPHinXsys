#ifndef UPDATE_CELL_LINKED_LIST_HPP
#define UPDATE_CELL_LINKED_LIST_HPP

#include "update_cell_linked_list.h"

#include "adaptation.hpp"
#include "base_particles.hpp"
#include "mesh_iterators.h"

namespace SPH
{
//=================================================================================================//
template <class MeshType>
ParticleCellLinkedList<MeshType>::ParticleCellLinkedList(
    const MeshType &mesh, Vecd *pos, UnsignedInt *particle_id_list,
    UnsignedInt *particle_offset_list)
    : MeshType(mesh), pos_(pos), particle_id_list_(particle_id_list),
      particle_offset_list_(particle_offset_list) {}
//=================================================================================================//
template <class MeshType>
template <typename NeighborhoodType, typename FunctionOnEach, typename... Args>
void ParticleCellLinkedList<MeshType>::forEachNeighbor(
    UnsignedInt index_i, const Vecd *source_pos,
    const NeighborhoodType &neighborhood, const FunctionOnEach &function, Args &&...args) const
{
    const Vecd pos_i = source_pos[index_i];
    const Arrayi target_cell_index = this->CellIndexFromPosition(pos_i);
    const int search_depth = neighborhood.get_search_depth(index_i);
    mesh_for_each(
        Arrayi::Zero().max(target_cell_index - search_depth * Arrayi::Ones()),
        this->all_cells_.min(target_cell_index + (search_depth + 1) * Arrayi::Ones()),
        [&](const Arrayi &cell_index)
        {
            const UnsignedInt linear_index = this->transferCellIndexTo1D(cell_index, this->all_cells_);
            // Since offset_cell_size_ has linear_cell_size_+1 elements, no boundary checks are needed.
            // offset_cell_size_[0] == 0 && offset_cell_size_[linear_cell_size_] == total_real_particles_
            for (UnsignedInt n = particle_offset_list_[linear_index]; n < particle_offset_list_[linear_index + 1]; ++n)
            {
                const UnsignedInt index_j = particle_id_list_[n];
                if (neighborhood.isWithinSupport(pos_i, pos_[index_j]))
                {
                    function(index_j, std::forward<Args>(args)...);
                }
            }
        });
}
//=================================================================================================//
template <typename MeshType>
template <class ExecutionPolicy>
UpdateCellLinkedList<MeshType>::
    UpdateCellLinkedList(const ExecutionPolicy &execution_policy, RealBody &real_body)
    : LocalDynamics(real_body),
      mesh_(createConstantEntity<MeshType>(
          execution_policy, "Mesh",
          real_body.sph_adaptation_->createBackGroundMesh<MeshType>(real_body))),
      number_of_cells_plus_one_(mesh_->NumberOfCells() + 1),
      particle_id_list_size_(SMAX(particles_->ParticlesBound(), number_of_cells_plus_one_)),
      pos_(particles_->getVariableDataByName<Vecd>(execution_policy, "Position")),
      particle_id_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleIDList", particle_id_list_size_)),
      particle_offset_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "ParticleOffsetList", number_of_cells_plus_one_)),
      current_size_list_(particles_->registerDiscreteVariable<UnsignedInt>(execution_policy, "CurrentCellSize", number_of_cells_plus_one_))
{
    particles_->addVariableToWrite<UnsignedInt>("ParticleIDList");
}
//=================================================================================================//
template <typename MeshType>
ParticleCellLinkedList<MeshType> UpdateCellLinkedList<MeshType>::getParticleCellLinkedList() const
{
    return ParticleCellLinkedList<MeshType>(mesh_, pos_, particle_id_list_, particle_offset_list_);
}
//=================================================================================================//
template <class MeshType, class ExecutionPolicy>
UpdateCellLinkedList<MeshType, ExecutionPolicy>::UpdateCellLinkedList(RealBody &real_body)
    : UpdateCellLinkedList<MeshType>(ExecutionPolicy{}, real_body),
      BaseDynamics<void>(){};
//=================================================================================================//
template <class MeshType, class ExecutionPolicy>
void UpdateCellLinkedList<MeshType, ExecutionPolicy>::exec(Real dt)
{
    this->clearParticleOffsetList(ExecutionPolicy{});
    this->incrementCellSize(ExecutionPolicy{});
    this->exclusiveScanParticleOffsetList(ExecutionPolicy{});
    this->updateCellLists(ExecutionPolicy{});
}
//=================================================================================================//
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_HPP
