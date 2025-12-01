#include "cell_linked_list.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void BaseCellLinkedList::tagBoundingCellsByMesh(Mesh &mesh, StdVec<CellLists> &cell_data_lists,
                                                const BoundingBoxd &bounding_bounds, int axis)
{
    int second_axis = NextAxis(axis);
    Array2i body_lower_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.lower_);
    Array2i body_upper_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.upper_);
    Array2i all_cells = mesh.AllCells();
    // lower bound cells
    for (int j = SMAX(body_lower_bound_cell_[second_axis] - 1, 0);
         j < SMIN(body_upper_bound_cell_[second_axis] + 2, all_cells[second_axis]); ++j)
        for (int i = SMAX(body_lower_bound_cell_[axis] - 1, 0);
             i < SMIN(body_lower_bound_cell_[axis] + 2, all_cells[axis]); ++i)
        {
            Array2i cell = Array2i::Zero();
            cell[axis] = i;
            cell[second_axis] = j;
            UnsignedInt linear_index = mesh.LinearCellIndex(cell);
            cell_data_lists[0].first.push_back(&cell_index_lists_[linear_index]);
            cell_data_lists[0].second.push_back(&cell_data_lists_[linear_index]);
        }

    // upper bound cells
    for (int j = SMAX(body_lower_bound_cell_[second_axis] - 1, 0);
         j < SMIN(body_upper_bound_cell_[second_axis] + 2, all_cells[second_axis]); ++j)
        for (int i = SMAX(body_upper_bound_cell_[axis] - 1, 0);
             i < SMIN(body_upper_bound_cell_[axis] + 2, all_cells[axis]); ++i)
        {
            Array2i cell = Array2i::Zero();
            cell[axis] = i;
            cell[second_axis] = j;
            UnsignedInt linear_index = mesh.LinearCellIndex(cell);
            cell_data_lists[1].first.push_back(&cell_index_lists_[linear_index]);
            cell_data_lists[1].second.push_back(&cell_data_lists_[linear_index]);
        }
}
//=================================================================================================//
} // namespace SPH
