#include "cell_linked_list.h"

#include "base_particles.hpp"
#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
void CellLinkedList ::allocateMeshDataMatrix()
{
    Allocate2dArray(cell_index_lists_, all_cells_);
    Allocate2dArray(cell_data_lists_, all_cells_);

    mesh_parallel_for(MeshRange(Array2i::Zero(), all_cells_),
                      [&](int i, int j)
                      {
                          cell_index_lists_[i][j].reserve(12);
                          cell_data_lists_[i][j].reserve(12);
                      });
}
//=================================================================================================//
void CellLinkedList ::deleteMeshDataMatrix()
{
    Delete2dArray(cell_index_lists_, all_cells_);
    Delete2dArray(cell_data_lists_, all_cells_);
}
//=================================================================================================//
void CellLinkedList::clearCellLists()
{
    mesh_parallel_for(MeshRange(Array2i::Zero(), all_cells_),
                      [&](int i, int j)
                      {
                          cell_index_lists_[i][j].clear();
                      });
}
//=================================================================================================//
void CellLinkedList::UpdateCellListData(BaseParticles &base_particles)
{
    StdLargeVec<Vecd> &pos = base_particles.pos_;
    StdLargeVec<Real> &Vol = base_particles.Vol_;
    mesh_parallel_for(
        MeshRange(Array2i::Zero(), all_cells_),
        [&](int i, int j)
        {
            cell_data_lists_[i][j].clear();
            ConcurrentIndexVector &cell_list = cell_index_lists_[i][j];
            for (size_t s = 0; s != cell_list.size(); ++s)
            {
                size_t index = cell_list[s];
                cell_data_lists_[i][j].emplace_back(std::make_tuple(index, pos[index], Vol[index]));
            }
        });
}
//=================================================================================================//
void CellLinkedList::updateSplitCellLists(SplitCellLists &split_cell_lists)
{
    // clear the data
    clearSplitCellLists(split_cell_lists);
    mesh_parallel_for(
        MeshRange(Array2i::Zero(), all_cells_),
        [&](int i, int j)
        {
            size_t real_particles_in_cell = cell_index_lists_[i][j].size();
            if (real_particles_in_cell != 0)
            {
                split_cell_lists[transferMeshIndexTo1D(Array2i(3, 3), Array2i(i % 3, j % 3))]
                    .push_back(&cell_index_lists_[i][j]);
            }
        });
}
//=================================================================================================//
void CellLinkedList ::insertParticleIndex(size_t particle_index, const Vecd &particle_position)
{
    Array2i cellpos = CellIndexFromPosition(particle_position);
    cell_index_lists_[cellpos[0]][cellpos[1]].emplace_back(particle_index);
}
//=================================================================================================//
void CellLinkedList ::InsertListDataEntry(
    size_t particle_index, const Vecd &particle_position, Real volumetric)
{
    Array2i cellpos = CellIndexFromPosition(particle_position);
    cell_data_lists_[cellpos[0]][cellpos[1]].emplace_back(
        std::make_tuple(particle_index, particle_position, volumetric));
}
//=================================================================================================//
ListData CellLinkedList::findNearestListDataEntry(const Vecd &position)
{
    Real min_distance_sqr = Infinity;
    ListData nearest_entry(MaxSize_t, Infinity * Vecd::Ones(), Infinity);

    Array2i cell = CellIndexFromPosition(position);
    mesh_for_each(
        Array2i::Zero().max(cell - Array2i::Ones()),
        all_cells_.min(cell + 2 * Array2i::Ones()),
        [&](int l, int m)
        {
            ListDataVector &target_particles = cell_data_lists_[l][m];
            for (const ListData &list_data : target_particles)
            {
                Real distance_sqr = (position - std::get<1>(list_data)).squaredNorm();
                if (distance_sqr < min_distance_sqr)
                {
                    min_distance_sqr = distance_sqr;
                    nearest_entry = list_data;
                }
            }
        });
    return nearest_entry;
}
//=================================================================================================//
void CellLinkedList::
    tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included)
{
    mesh_parallel_for(
        MeshRange(Array2i::Zero(), all_cells_),
        [&](int i, int j)
        {
            bool is_included = false;
            mesh_for_each(
                Array2i::Zero().max(Array2i(i, j) - Array2i::Ones()),
                all_cells_.min(Array2i(i, j) + 2 * Array2i::Ones()),
                [&](int l, int m)
                {
                    if (check_included(CellPositionFromIndex(Array2i(l, m)), grid_spacing_))
                    {
                        is_included = true;
                    }
                });
            if (is_included == true)
                cell_lists.push_back(&cell_index_lists_[i][j]);
        });
}
//=================================================================================================//
void CellLinkedList::
    tagBoundingCells(StdVec<CellLists> &cell_data_lists, BoundingBox &bounding_bounds, int axis)
{
    int second_axis = SecondAxis(axis);
    Array2i body_lower_bound_cell_ = CellIndexFromPosition(bounding_bounds.first_);
    Array2i body_upper_bound_cell_ = CellIndexFromPosition(bounding_bounds.second_);

    // lower bound cells
    for (int j = SMAX(body_lower_bound_cell_[second_axis] - 1, 0);
         j < SMIN(body_upper_bound_cell_[second_axis] + 2, all_cells_[second_axis]); ++j)
        for (int i = SMAX(body_lower_bound_cell_[axis] - 1, 0);
             i < SMIN(body_lower_bound_cell_[axis] + 2, all_cells_[axis]); ++i)
        {
            Array2i cell = Array2i::Zero();
            cell[axis] = i;
            cell[second_axis] = j;
            cell_data_lists[0].first.push_back(&cell_index_lists_[cell[0]][cell[1]]);
            cell_data_lists[0].second.push_back(&cell_data_lists_[cell[0]][cell[1]]);
        }

    // upper bound cells
    for (int j = SMAX(body_lower_bound_cell_[second_axis] - 1, 0);
         j < SMIN(body_upper_bound_cell_[second_axis] + 2, all_cells_[second_axis]); ++j)
        for (int i = SMAX(body_upper_bound_cell_[axis] - 1, 0);
             i < SMIN(body_upper_bound_cell_[axis] + 2, all_cells_[axis]); ++i)
        {
            Array2i cell = Array2i::Zero();
            cell[axis] = i;
            cell[second_axis] = j;
            cell_data_lists[1].first.push_back(&cell_index_lists_[cell[0]][cell[1]]);
            cell_data_lists[1].second.push_back(&cell_data_lists_[cell[0]][cell[1]]);
        }
}
//=============================================================================================//
void CellLinkedList::writeMeshFieldToPlt(std::ofstream &output_file)
{
    Array2i number_of_operation = all_cells_;

    output_file << "\n";
    output_file << "title='View'"
                << "\n";
    output_file << "variables= "
                << "x, "
                << "y, "
                << "particles_in_cell "
                << "\n";
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
                << "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            Vecd data_position = CellPositionFromIndex(Array2i(i, j));
            output_file << data_position[0] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            Vecd data_position = CellPositionFromIndex(Array2i(i, j));
            output_file << data_position[1] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            output_file << cell_index_lists_[i][j].size() << " ";
        }
        output_file << " \n";
    }
}
//=================================================================================================//
} // namespace SPH
