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
void CellLinkedList::particle_for_split(const execution::SequencedPolicy &, const LocalDynamicsFunction &local_dynamics_function)
{
    // foward sweeping
    for (size_t k = 0; k < number_of_split_cell_lists_; k++)
    {
        const Arrayi split_cell_index = transfer1DtoMeshIndex(3 * Arrayi::Ones(), k);
        const Arrayi all_cells_k = (all_cells_ - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const size_t number_of_cells = get1DMeshSize(all_cells_k);

        for (size_t l = 0; l < number_of_cells; l++)
        {
            const Arrayi cell_index = split_cell_index + 3 * transfer1DtoMeshIndex(all_cells_k, l);
            const ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
            for (const size_t index_i : cell_list)
            {
                local_dynamics_function(index_i);
            }
        }
    }

    // backward sweeping
    for (size_t k = number_of_split_cell_lists_; k != 0; --k)
    {
        const Arrayi split_cell_index = transfer1DtoMeshIndex(3 * Arrayi::Ones(), k - 1);
        const Arrayi all_cells_k = (all_cells_ - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const size_t number_of_cells = get1DMeshSize(all_cells_k);

        for (size_t l = 0; l < number_of_cells; l++)
        {
            const Arrayi cell_index = split_cell_index + 3 * transfer1DtoMeshIndex(all_cells_k, l);
            const ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
            for (size_t i = cell_list.size(); i != 0; --i)
            {
                local_dynamics_function(cell_list[i - 1]);
            }
        }
    }
}
//=================================================================================================//
template <class LocalDynamicsFunction>
void CellLinkedList::particle_for_split(const execution::ParallelPolicy &, const LocalDynamicsFunction &local_dynamics_function)
{
    // foward sweeping
    for (size_t k = 0; k < number_of_split_cell_lists_; k++)
    {
        const Arrayi split_cell_index = transfer1DtoMeshIndex(3 * Arrayi::Ones(), k);
        const Arrayi all_cells_k = (all_cells_ - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const size_t number_of_cells = get1DMeshSize(all_cells_k);

        parallel_for(
            IndexRange(0, number_of_cells),
            [&](const IndexRange &r)
            {
                for (size_t l = r.begin(); l < r.end(); ++l)
                {
                    const Arrayi cell_index = split_cell_index + 3 * transfer1DtoMeshIndex(all_cells_k, l);
                    const ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
                    for (const size_t index_i : cell_list)
                    {
                        local_dynamics_function(index_i);
                    }
                }
            },
            ap);
    }

    // backward sweeping
    for (size_t k = number_of_split_cell_lists_; k != 0; --k)
    {
        const Arrayi split_cell_index = transfer1DtoMeshIndex(3 * Arrayi::Ones(), k - 1);
        const Arrayi all_cells_k = (all_cells_ - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const size_t number_of_cells = get1DMeshSize(all_cells_k);

        parallel_for(
            IndexRange(0, number_of_cells),
            [&](const IndexRange &r)
            {
                for (size_t l = r.begin(); l < r.end(); ++l)
                {
                    const Arrayi cell_index = split_cell_index + 3 * transfer1DtoMeshIndex(all_cells_k, l);
                    const ConcurrentIndexVector &cell_list = getCellDataList(cell_index_lists_, cell_index);
                    for (size_t i = cell_list.size(); i != 0; --i)
                    {
                        local_dynamics_function(cell_list[i - 1]);
                    }
                }
            },
            ap);
    }
}
//=================================================================================================//
} // namespace SPH
