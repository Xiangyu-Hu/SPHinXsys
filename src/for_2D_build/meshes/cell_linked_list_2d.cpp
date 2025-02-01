#include "cell_linked_list.h"

#include "base_particles.hpp"

namespace SPH
{
//=============================================================================================//
void BaseCellLinkedList::writeMeshFieldToPltByMesh(Mesh &mesh, size_t mesh_offset, std::ofstream &output_file)
{
    Array2i number_of_operation = mesh.AllCells();

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
            Vecd data_position = mesh.CellPositionFromIndex(Array2i(i, j));
            output_file << data_position[0] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            Vecd data_position = mesh.CellPositionFromIndex(Array2i(i, j));
            output_file << data_position[1] << " ";
        }
        output_file << " \n";
    }

    for (int j = 0; j != number_of_operation[1]; ++j)
    {
        for (int i = 0; i != number_of_operation[0]; ++i)
        {
            size_t linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(Array2i(i, j));
            output_file << cell_index_lists_[linear_index].size() << " ";
        }
        output_file << " \n";
    }
}
//=================================================================================================//
void BaseCellLinkedList::tagBoundingCellsByMesh(Mesh &mesh, size_t mesh_offset,
                                                StdVec<CellLists> &cell_data_lists,
                                                const BoundingBox &bounding_bounds, int axis)
{
    int second_axis = NextAxis(axis);
    Array2i body_lower_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.first_);
    Array2i body_upper_bound_cell_ = mesh.CellIndexFromPosition(bounding_bounds.second_);
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
            size_t linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell);
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
            size_t linear_index = mesh_offset + mesh.LinearCellIndexFromCellIndex(cell);
            cell_data_lists[1].first.push_back(&cell_index_lists_[linear_index]);
            cell_data_lists[1].second.push_back(&cell_data_lists_[linear_index]);
        }
}
//=================================================================================================//
} // namespace SPH
