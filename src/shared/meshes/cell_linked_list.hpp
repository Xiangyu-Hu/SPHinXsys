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
void CellLinkedList::searchNeighborsByParticles(
    DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
    GetSearchDepth &get_search_depth, GetNeighborRelation &get_neighbor_relation)
{
    StdLargeVec<Vecd> &pos = dynamics_range.getBaseParticles().ParticlePositions();
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
template <class LocalDynamicsFunction>
void CellLinkedList::particle_for_split(const LocalDynamicsFunction &local_dynamics_function)
{
    mesh_split_parallel_for(
        MeshRange(Arrayi::Zero(), all_cells_),
        3 * Arrayi::Ones(),
        [&](const Arrayi &cell_index)
        {
            const ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
            for (size_t index_i : cell_list)
            {
                local_dynamics_function(index_i);
            }
        });
}
//=================================================================================================//
} // namespace SPH
