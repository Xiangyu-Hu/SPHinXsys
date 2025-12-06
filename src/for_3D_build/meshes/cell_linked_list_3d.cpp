/**
 * @file 	cell_linked_list_supplementary.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "cell_linked_list.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void BaseCellLinkedList::tagBoundingCellsByMesh(Mesh &mesh, StdVec<CellLists> &cell_data_lists,
                                                const BoundingBoxd &bounding_bounds, int axis)
{
    int second_axis = NextAxis(axis);
    int third_axis = NextNextAxis(axis);
    Array3i body_lower_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.lower_);
    Array3i body_upper_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.upper_);
    Array3i all_cells = mesh.AllCells();

    // lower bound cells
    for (int k = SMAX(body_lower_bound_cell_[third_axis] - 1, 0);
         k < SMIN(body_upper_bound_cell_[third_axis] + 2, all_cells[third_axis]); ++k)
    {
        for (int j = SMAX(body_lower_bound_cell_[second_axis] - 1, 0);
             j < SMIN(body_upper_bound_cell_[second_axis] + 2, all_cells[second_axis]); ++j)
        {
            for (int i = SMAX(body_lower_bound_cell_[axis] - 1, 0);
                 i < SMIN(body_lower_bound_cell_[axis] + 2, all_cells[axis]); ++i)
            {
                Array3i cell = Array3i::Zero();
                cell[axis] = i;
                cell[second_axis] = j;
                cell[third_axis] = k;
                UnsignedInt linear_index = mesh.LinearCellIndex(cell);
                cell_data_lists[0].first.push_back(&cell_index_lists_[linear_index]);
                cell_data_lists[0].second.push_back(&cell_data_lists_[linear_index]);
            }
        }
    }

    // upper bound cells
    for (int k = SMAX(body_lower_bound_cell_[third_axis] - 1, 0);
         k < SMIN(body_upper_bound_cell_[third_axis] + 2, all_cells[third_axis]); ++k)
    {
        for (int j = SMAX(body_lower_bound_cell_[second_axis] - 1, 0);
             j < SMIN(body_upper_bound_cell_[second_axis] + 2, all_cells[second_axis]); ++j)
        {
            for (int i = SMAX(body_upper_bound_cell_[axis] - 1, 0);
                 i < SMIN(body_upper_bound_cell_[axis] + 2, all_cells[axis]); ++i)
            {
                Array3i cell = Array3i::Zero();
                cell[axis] = i;
                cell[second_axis] = j;
                cell[third_axis] = k;
                UnsignedInt linear_index = mesh.LinearCellIndex(cell);
                cell_data_lists[1].first.push_back(&cell_index_lists_[linear_index]);
                cell_data_lists[1].second.push_back(&cell_data_lists_[linear_index]);
            }
        }
    }
}
//=================================================================================================//
} // namespace SPH
