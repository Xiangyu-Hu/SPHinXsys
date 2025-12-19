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
template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
void BaseCellLinkedList::searchNeighborsByMesh(
    Mesh &mesh, DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
    GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation)
{
    Vecd *pos = dynamics_range.getBaseParticles().ParticlePositions();
    particle_for(execution::ParallelPolicy(), dynamics_range.LoopRange(),
                 [&](UnsignedInt index_i)
                 {
                     int search_depth = get_search_depth(index_i);
                     Arrayi target_cell_index = mesh.CellIndexFromPosition(pos[index_i]);

                     Neighborhood &neighborhood = particle_configuration[index_i];
                     mesh_for_each(
                         Arrayi::Zero().max(target_cell_index - search_depth * Arrayi::Ones()),
                         mesh.AllCells().min(target_cell_index + (search_depth + 1) * Arrayi::Ones()),
                         [&](const Arrayi &cell_index)
                         {
                             UnsignedInt linear_index = mesh.LinearCellIndex(cell_index);
                             ListDataVector &target_particles = cell_data_lists_[linear_index];
                             for (const ListData &data_list : target_particles)
                             {
                                 get_neighbor_relation(neighborhood, pos[index_i], index_i, data_list);
                             }
                         });
                 });
}
//=================================================================================================//
template <class ExecutionPolicy, class LocalDynamicsFunction>
void BaseCellLinkedList::particle_for_split_by_mesh(
    const ExecutionPolicy &ex_policy, Mesh &mesh, const LocalDynamicsFunction &local_dynamics_function)
{
    const Arrayi array3s = 3 * Arrayi::Ones();

    // forward sweeping
    for (int k = 0; k < array3s.prod(); k++)
    {
        // get the corresponding 2D/3D split cell index (m, n)
        // e.g., for k = 0, split_cell_index = (0,0), for k = 3, split_cell_index = (1,0), etc.
        const Arrayi split_cell_index = Mesh::transfer1DtoMeshIndex(array3s, k);
        // get the number of cells belonging to the split cell k
        // i_max = (M - m - 1) / 3 + 1, j_max = (N - n - 1) / 3 + 1
        // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), then i_max = 2, j_max = 3
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();

        // looping over all cells in the split cell k
        particle_for(ex_policy, IndexRange(0, all_cells_k.prod()),
                     [&](UnsignedInt l)
                     {
                         // get the 2D/3D cell index of the l-th cell in the split cell k
                         // (i , j) = (m + 3 * (l / j_max), n + 3 * l % i_max)
                         // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), l = 0, then (i, j) = (1, 1)
                         // l = 1, then (i, j) = (1, 4), l = 3, then (i, j) = (4, 1), etc.
                         const Arrayi cell_index = split_cell_index + 3 * Mesh::transfer1DtoMeshIndex(all_cells_k, l);
                         UnsignedInt linear_index = mesh.LinearCellIndex(cell_index);
                         // get the list of particles in the cell (i, j)
                         const ConcurrentIndexVector &cell_list = cell_index_lists_[linear_index];
                         // looping over all particles in the cell (i, j)
                         for (const UnsignedInt index_i : cell_list)
                         {
                             local_dynamics_function(index_i);
                         }
                     });
    }

    // backward sweeping
    for (int k = array3s.prod(); k != 0; --k)
    {
        const Arrayi split_cell_index = Mesh::transfer1DtoMeshIndex(array3s, k - 1);
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();

        particle_for(ex_policy, IndexRange(0, all_cells_k.prod()),
                     [&](UnsignedInt l)
                     {
                         const Arrayi cell_index = split_cell_index + 3 * Mesh::transfer1DtoMeshIndex(all_cells_k, l);
                         UnsignedInt linear_index = mesh.LinearCellIndex(cell_index);
                         const ConcurrentIndexVector &cell_list = cell_index_lists_[linear_index];
                         for (UnsignedInt i = cell_list.size(); i != 0; --i)
                         {
                             local_dynamics_function(cell_list[i - 1]);
                         }
                     });
    }
}
//=================================================================================================//
template <class ExecutionPolicy, class LocalDynamicsFunction>
void BaseCellLinkedList::particle_for_split(const ExecutionPolicy &ex_policy,
                                            const LocalDynamicsFunction &local_dynamics_function)
{
    for (UnsignedInt level = 0; level != resolution_levels_; ++level)
        particle_for_split_by_mesh(ex_policy, getMesh(level), local_dynamics_function);
}
//=================================================================================================//
template <class ExecutionPolicy>
NeighborSearch::NeighborSearch(const ExecutionPolicy &ex_policy, CellLinkedList &cell_linked_list)
    : Mesh(cell_linked_list.getMesh()),
      particle_index_(cell_linked_list.dvParticleIndex()->DelegatedData(ex_policy)),
      cell_offset_(cell_linked_list.dvCellOffset()->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename FunctionOnEach>
void NeighborSearch::forEachSearch(const Vecd &source_pos, const FunctionOnEach &function,
                                   const BoundingBoxi &search_box) const
{
    const BoundingBoxi search_range =
        search_box.translate(CellIndexFromPosition(source_pos));
    mesh_for_each(
        Arrayi::Zero().max(search_range.lower_), all_cells_.min(search_range.upper_ + Arrayi::Ones()),
        [&](const Arrayi &cell_index)
        {
            const UnsignedInt linear_index = LinearCellIndex(cell_index);
            // Since offset_cell_size_ has linear_cell_size_+1 elements, no boundary checks are needed.
            // offset_cell_size_[0] == 0 && offset_cell_size_[linear_cell_size_] == total_real_particles_
            for (UnsignedInt n = cell_offset_[linear_index]; n < cell_offset_[linear_index + 1]; ++n)
            {
                function(particle_index_[n]);
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
} // namespace SPH
