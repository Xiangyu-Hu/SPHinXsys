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
    Mesh &mesh, UnsignedInt mesh_offset,
    DynamicsRange &dynamics_range, ParticleConfiguration &particle_configuration,
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
                             UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
                             ListDataVector &target_particles = cell_data_lists_[linear_index];
                             for (const ListData &data_list : target_particles)
                             {
                                 get_neighbor_relation(neighborhood, pos[index_i], index_i, data_list);
                             }
                         });
                 });
}
//=================================================================================================//
template <class LocalDynamicsFunction>
void BaseCellLinkedList::particle_for_split_by_mesh(
    const execution::SequencedPolicy &, Mesh &mesh, UnsignedInt mesh_offset,
    const LocalDynamicsFunction &local_dynamics_function)
{
    // forward sweeping
    for (UnsignedInt k = 0; k < number_of_split_cell_lists_; k++)
    {
        // get the corresponding 2D/3D split cell index (m, n)
        // e.g., for k = 0, split_cell_index = (0,0), for k = 3, split_cell_index = (1,0), etc.
        const Arrayi split_cell_index = mesh.transfer1DtoMeshIndex(3 * Arrayi::Ones(), k);
        // get the number of cells belonging to the split cell k
        // i_max = (M - m - 1) / 3 + 1, j_max = (N - n - 1) / 3 + 1
        // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), then i_max = 2, j_max = 3
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const UnsignedInt number_of_cells = all_cells_k.prod(); // i_max * j_max

        // looping over all cells in the split cell k
        for (UnsignedInt l = 0; l < number_of_cells; l++)
        {
            // get the 2D/3D cell index of the l-th cell in the split cell k
            // (i , j) = (m + 3 * (l / j_max), n + 3 * l % i_max)
            // e.g. all_cells = (M,N) = (6, 9), (m, n) = (1, 1), l = 0, then (i, j) = (1, 1)
            // l = 1, then (i, j) = (1, 4), l = 3, then (i, j) = (4, 1), etc.
            const Arrayi cell_index = split_cell_index + 3 * mesh.transfer1DtoMeshIndex(all_cells_k, l);
            UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
            // get the list of particles in the cell (i, j)
            const ConcurrentIndexVector &cell_list = cell_index_lists_[linear_index];
            // looping over all particles in the cell (i, j)
            for (const UnsignedInt index_i : cell_list)
            {
                local_dynamics_function(index_i);
            }
        }
    }

    // backward sweeping
    for (UnsignedInt k = number_of_split_cell_lists_; k != 0; --k)
    {
        const Arrayi split_cell_index = mesh.transfer1DtoMeshIndex(3 * Arrayi::Ones(), k - 1);
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const UnsignedInt number_of_cells = all_cells_k.prod();

        for (UnsignedInt l = 0; l < number_of_cells; l++)
        {
            const Arrayi cell_index = split_cell_index + 3 * mesh.transfer1DtoMeshIndex(all_cells_k, l);
            UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
            const ConcurrentIndexVector &cell_list = cell_index_lists_[linear_index];
            for (UnsignedInt i = cell_list.size(); i != 0; --i)
            {
                local_dynamics_function(cell_list[i - 1]);
            }
        }
    }
}
//=================================================================================================//
template <class LocalDynamicsFunction>
void BaseCellLinkedList::particle_for_split_by_mesh(
    const execution::ParallelPolicy &, Mesh &mesh, UnsignedInt mesh_offset,
    const LocalDynamicsFunction &local_dynamics_function)
{
    // forward sweeping
    for (UnsignedInt k = 0; k < number_of_split_cell_lists_; k++)
    {
        const Arrayi split_cell_index = mesh.transfer1DtoMeshIndex(3 * Arrayi::Ones(), k);
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const UnsignedInt number_of_cells = all_cells_k.prod();

        parallel_for(
            IndexRange(0, number_of_cells),
            [&](const IndexRange &r)
            {
                for (UnsignedInt l = r.begin(); l < r.end(); ++l)
                {
                    const Arrayi cell_index = split_cell_index + 3 * mesh.transfer1DtoMeshIndex(all_cells_k, l);
                    UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
                    const ConcurrentIndexVector &cell_list = cell_index_lists_[linear_index];
                    for (const UnsignedInt index_i : cell_list)
                    {
                        local_dynamics_function(index_i);
                    }
                }
            },
            ap);
    }

    // backward sweeping
    for (UnsignedInt k = number_of_split_cell_lists_; k != 0; --k)
    {
        const Arrayi split_cell_index = mesh.transfer1DtoMeshIndex(3 * Arrayi::Ones(), k - 1);
        const Arrayi all_cells_k = (mesh.AllCells() - split_cell_index - Arrayi::Ones()) / 3 + Arrayi::Ones();
        const UnsignedInt number_of_cells = all_cells_k.prod();

        parallel_for(
            IndexRange(0, number_of_cells),
            [&](const IndexRange &r)
            {
                for (UnsignedInt l = r.begin(); l < r.end(); ++l)
                {
                    const Arrayi cell_index = split_cell_index + 3 * mesh.transfer1DtoMeshIndex(all_cells_k, l);
                    UnsignedInt linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell_index);
                    const ConcurrentIndexVector &cell_list = cell_index_lists_[linear_index];
                    for (UnsignedInt i = cell_list.size(); i != 0; --i)
                    {
                        local_dynamics_function(cell_list[i - 1]);
                    }
                }
            },
            ap);
    }
}
//=================================================================================================//
template <typename DataType>
DataType *BaseCellLinkedList::initializeVariable(DiscreteVariable<DataType> *variable, DataType initial_value)
{
    DataType *data_field = variable->Data();
    for (size_t i = 0; i != variable->getDataSize(); ++i)
    {
        data_field[i] = initial_value;
    }
    return data_field;
}
//=================================================================================================//
template <typename DataType, typename... Args>
DiscreteVariable<DataType> *BaseCellLinkedList::registerDiscreteVariableOnly(
    const std::string &name, size_t data_size, Args &&...args)
{
    DiscreteVariable<DataType> *variable = findVariableByName<DataType>(all_discrete_variables_, name);
    if (variable == nullptr)
    {
        variable = addVariableToAssemble<DataType>(all_discrete_variables_, all_discrete_variable_ptrs_,
                                                   name, data_size);
        initializeVariable(variable, std::forward<Args>(args)...);
    }
    return variable;
}
//=================================================================================================//
template <class ExecutionPolicy>
NeighborSearch::NeighborSearch(const ExecutionPolicy &ex_policy, CellLinkedList &cell_linked_list)
    : Mesh(cell_linked_list.getMesh()),
      particle_index_(cell_linked_list.dvParticleIndex()->DelegatedData(ex_policy)),
      cell_offset_(cell_linked_list.dvCellOffset()->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename FunctionOnEach>
void NeighborSearch::forEachSearch(UnsignedInt source_index, const Vecd *source_pos,
                                   const FunctionOnEach &function) const
{
    const Arrayi target_cell_index = CellIndexFromPosition(source_pos[source_index]);
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
template <class LocalDynamicsFunction>
void CellLinkedList::particle_for_split(const execution::SequencedPolicy &,
                                        const LocalDynamicsFunction &local_dynamics_function)
{
    particle_for_split_by_mesh(execution::SequencedPolicy(), *mesh_, 0, local_dynamics_function);
}
//=================================================================================================//
template <class LocalDynamicsFunction>
void CellLinkedList::particle_for_split(const execution::ParallelPolicy &,
                                        const LocalDynamicsFunction &local_dynamics_function)
{
    particle_for_split_by_mesh(execution::ParallelPolicy(), *mesh_, 0, local_dynamics_function);
}
//=================================================================================================//
template <class LocalDynamicsFunction>
void MultilevelCellLinkedList::particle_for_split(const execution::SequencedPolicy &seq,
                                                  const LocalDynamicsFunction &local_dynamics_function)
{
    for (UnsignedInt level = 0; level != meshes_.size(); ++level)
        particle_for_split_by_mesh(execution::SequencedPolicy(),
                                   *meshes_[level], mesh_offsets_[level], local_dynamics_function);
}
//=================================================================================================//
template <class LocalDynamicsFunction>
void MultilevelCellLinkedList::particle_for_split(const execution::ParallelPolicy &par,
                                                  const LocalDynamicsFunction &local_dynamics_function)
{
    for (UnsignedInt level = 0; level != meshes_.size(); ++level)
        particle_for_split_by_mesh(execution::ParallelPolicy(),
                                   *meshes_[level], mesh_offsets_[level], local_dynamics_function);
}
//=================================================================================================//
} // namespace SPH
