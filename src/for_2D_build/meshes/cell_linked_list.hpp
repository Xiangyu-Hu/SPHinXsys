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
    StdLargeVec<Vecd> &pos = dynamics_range.getBaseParticles().pos_;
    particle_for(execution::ParallelPolicy(), dynamics_range.LoopRange(),
                 [&](size_t index_i)
                 {
                     int search_depth = get_search_depth(index_i);
                     Array2i target_cell_index = CellIndexFromPosition(pos[index_i]);

                     Neighborhood &neighborhood = particle_configuration[index_i];
                     mesh_for_each(
                         Array2i::Zero().max(target_cell_index - search_depth * Array2i::Ones()),
                         all_cells_.min(target_cell_index + (search_depth + 1) * Array2i::Ones()),
                         [&](int l, int m)
                         {
                             ListDataVector &target_particles = cell_data_lists_[l][m];
                             for (const ListData &list_data : target_particles)
                             {
                                 get_neighbor_relation(neighborhood, pos[index_i], index_i, list_data);
                             }
                         });
                 });
}
//=================================================================================================//
} // namespace SPH
