/**
 * @file 	cell_linked_list.hpp
 * @brief 	Here gives the classes for managing cell linked lists. This is the basic class
 * 			for building the particle configurations.
 * @details The cell linked list saves for each body a list of particles
 * 			located within the cell.
 * @author	Chi Zhang, Yongchuan and Xiangyu Hu
 */

#pragma once

#include "base_particles.h"
#include "cell_linked_list.h"
#include "mesh_iterators.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
NeighborSearch::NeighborSearch(const ExecutionPolicy &ex_policy, CellLinkedList &cell_linked_list)
    : Mesh(cell_linked_list), grid_spacing_squared_(grid_spacing_ * grid_spacing_),
      pos_(cell_linked_list.getParticlePosition()->DelegatedDataField(ex_policy)),
      particle_index_(cell_linked_list.getParticleIndex()->DelegatedDataField(ex_policy)),
      cell_offset_(cell_linked_list.getCellOffset()->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <typename FunctionOnEach>
void NeighborSearch::forEachSearch(UnsignedInt index_i, const Vecd *source_pos,
                                   const FunctionOnEach &function) const
{
    const Arrayi target_cell_index = CellIndexFromPosition(source_pos[index_i]);
    mesh_for_each(
        Arrayi::Zero().max(target_cell_index - Arrayi::Ones()),
        all_cells_.min(target_cell_index + 2 * Arrayi::Ones()),
        [&](const Arrayi &cell_index)
        {
            const UnsignedInt linear_index = LinearCellIndexFromCellIndex(cell_index);
            // Since offset_cell_size_ has linear_cell_size_+1 elements, no boundary checks are needed.
            // offset_cell_size_[0] == 0 && offset_cell_size_[linear_cell_size_] == total_real_particles_
            for (UnsignedInt n = cell_offset_[linear_index]; n < cell_offset_[linear_index + 1]; ++n)
            {
                const UnsignedInt index_j = particle_index_[n];
                if ((source_pos[index_i] - pos_[index_j]).squaredNorm() < grid_spacing_squared_)
                {
                    function(index_j);
                }
            }
        });
}
//=================================================================================================//
template <class ExecutionPolicy>
NeighborSearch CellLinkedList::createNeighborSearch(const ExecutionPolicy &ex_policy)
{
    return NeighborSearch(ex_policy, *this);
}
//=================================================================================================//
template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
void CellLinkedList::searchNeighborsByParticles(
    DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
    GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation)
{
    Vecd *pos = dynamics_range.getBaseParticles().ParticlePositions();
    particle_for(execution::ParallelPolicy(), dynamics_range.LoopRange(),
                 [&](size_t index_i)
                 {
                     int search_depth = get_search_depth(index_i);
                     Arrayi target_cell_index = CellIndexFromPosition(pos[index_i]);

                     Neighborhood &neighborhood = particle_configuration[index_i];
                     mesh_for_each(
                         Arrayi::Zero().max(target_cell_index - search_depth * Arrayi::Ones()),
                         all_cells_.min(target_cell_index + (search_depth + 1) * Arrayi::Ones()),
                         [&](const Arrayi &cell_index)
                         {
                             ListDataVector &target_particles = getCellDataList(cell_data_lists_, cell_index);
                             for (const ListData &data_list : target_particles)
                             {
                                 get_neighbor_relation(neighborhood, pos[index_i], index_i, data_list);
                             }
                         });
                 });
}
//=================================================================================================//
} // namespace SPH
