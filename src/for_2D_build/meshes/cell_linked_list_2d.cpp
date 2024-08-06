#include "cell_linked_list.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void CellLinkedList::
    tagBoundingCells(StdVec<CellLists> &cell_data_lists, const BoundingBox &bounding_bounds, int axis)
{
    int second_axis = NextAxis(axis);
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
            cell_data_lists[0].first.push_back(&getCellDataList(cell_index_lists_, cell));
            cell_data_lists[0].second.push_back(&getCellDataList(cell_data_lists_, cell));
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
            cell_data_lists[1].first.push_back(&getCellDataList(cell_index_lists_, cell));
            cell_data_lists[1].second.push_back(&getCellDataList(cell_data_lists_, cell));
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
            output_file << getCellDataList(cell_index_lists_, Array2i(i, j)).size() << " ";
        }
        output_file << " \n";
    }
}
//=================================================================================================//
} // namespace SPH
