/**
 * @file 	cell_linked_list_supplementary.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "cell_linked_list.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void BaseCellLinkedList::writeMeshFieldToPltByMesh(Mesh &mesh, std::ofstream &output_file)
{
    Array3i number_of_operation = mesh.AllCells();

    output_file << "\n";
    output_file << "title='View'"
                << "\n";
    output_file << "variables= "
                << "x, "
                << "y, "
                << "z, "
                << "particles_in_cell "
                << "\n";
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
                << "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = mesh.CellPositionFromIndex(Array3i(i, j, k));
                output_file << data_position[0] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = mesh.CellPositionFromIndex(Array3i(i, j, k));
                output_file << data_position[1] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = mesh.CellPositionFromIndex(Array3i(i, j, k));
                output_file << data_position[2] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                UnsignedInt linear_index = mesh.LinearCellIndex(Array3i(i, j, k));
                output_file << cell_data_lists_[linear_index].size() << " ";
            }
            output_file << " \n";
        }
}
//=================================================================================================//
void BaseCellLinkedList::tagBoundingCellsByMesh(Mesh &mesh, StdVec<CellLists> &cell_data_lists,
                                                const BoundingBox &bounding_bounds, int axis)
{
    int second_axis = NextAxis(axis);
    int third_axis = NextNextAxis(axis);
    Array3i body_lower_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.first_);
    Array3i body_upper_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.second_);
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
