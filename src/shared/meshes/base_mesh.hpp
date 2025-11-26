#ifndef BASE_MESH_HPP
#define BASE_MESH_HPP

#include "base_mesh.h"

namespace SPH
{
//=================================================================================================//
inline Arrayi Mesh::CellIndexFromPosition(const Vecd &position) const
{
    return floor((position - mesh_lower_bound_).array() / grid_spacing_)
        .cast<int>()
        .max(Arrayi::Zero())
        .min(all_grid_points_ - 2 * Arrayi::Ones());
}
//=================================================================================================//
inline UnsignedInt Mesh::LinearCellIndexFromPosition(const Vecd &position) const
{
    return linear_cell_index_offset_ +
           transferMeshIndexTo1D(all_cells_, CellIndexFromPosition(position));
}
//=================================================================================================//
inline UnsignedInt Mesh::LinearCellIndex(const Arrayi &cell_index) const
{
    return linear_cell_index_offset_ + transferMeshIndexTo1D(all_cells_, cell_index);
}
//=================================================================================================//
inline Vecd Mesh::CellLowerCornerPosition(const Arrayi &cell_index) const
{
    return mesh_lower_bound_ + cell_index.cast<Real>().matrix() * grid_spacing_;
}
//=================================================================================================//
inline Arrayi Mesh::boundCellIndex(const Arrayi &input) const
{
    Arrayi output = input;
    for (int i = 0; i < Dimensions; ++i)
    {
        if (output[i] < 0)
            output[i] = 0;
        if (output[i] >= all_cells_[i])
            output[i] = all_cells_[i] - 1;
    }
    return output;
}
//=================================================================================================//
inline Arrayi Mesh::DimensionalCellIndex(UnsignedInt linear_index) const
{
    return transfer1DtoMeshIndex(all_cells_, linear_index - linear_cell_index_offset_);
}
//=================================================================================================//
inline Array2i Mesh::transfer1DtoMeshIndex(const Array2i &mesh_size, UnsignedInt i)
{
    UnsignedInt row_size = mesh_size[1];
    UnsignedInt column = i / row_size;
    return Array2i(column, i - column * row_size);
}
//=================================================================================================//
inline Array3i Mesh::transfer1DtoMeshIndex(const Array3i &mesh_size, UnsignedInt i)
{
    UnsignedInt row_times_column_size = mesh_size[1] * mesh_size[2];
    UnsignedInt page = i / row_times_column_size;
    UnsignedInt left_over = (i - page * row_times_column_size);
    UnsignedInt row_size = mesh_size[2];
    UnsignedInt column = left_over / row_size;
    return Array3i(page, column, left_over - column * row_size);
}
//=================================================================================================//
inline UnsignedInt Mesh::transferMeshIndexTo1D(const Array2i &mesh_size, const Array2i &mesh_index)
{
    return mesh_index[0] * mesh_size[1] + mesh_index[1];
}
//=================================================================================================//
inline UnsignedInt Mesh::transferMeshIndexTo1D(const Array3i &mesh_size, const Array3i &mesh_index)
{
    return mesh_index[0] * mesh_size[1] * mesh_size[2] +
           mesh_index[1] * mesh_size[2] +
           mesh_index[2];
}
//=================================================================================================//
inline UnsignedInt Mesh::transferMeshIndexToMortonOrder(const Array2i &mesh_index)
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
}
//=================================================================================================//
inline UnsignedInt Mesh::transferMeshIndexToMortonOrder(const Array3i &mesh_index)
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1) | (MortonCode(mesh_index[2]) << 2);
}
//=================================================================================================//
inline UnsignedInt Mesh::MortonCode(const UnsignedInt &i)
{
    UnsignedInt x = i;
    x &= 0x3ff;
    x = (x | x << 16) & 0x30000ff;
    x = (x | x << 8) & 0x300f00f;
    x = (x | x << 4) & 0x30c30c3;
    x = (x | x << 2) & 0x9249249;
    return x;
}
//=============================================================================================//
template <typename DataType>
DiscreteVariable<DataType> *MultiLevelMeshField::getCellVariable(
    const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable =
        findVariableByName<DataType, DiscreteVariable>(all_cell_variables_, variable_name);
    if (variable == nullptr)
    {
        std::cout << "\n Error: the cell variable '" << variable_name << "' is not exist!" << std::endl;
        exit(1);
    }
    return variable;
}
//=============================================================================================//
template <typename DataType>
void MultiLevelMeshField::addCellVariableToWrite(const std::string &variable_name)
{
    DiscreteVariable<DataType> *variable = getCellVariable<DataType>(variable_name);
    addVariableToList<DiscreteVariable, DataType>(all_cell_variables_, variable);
}
//=============================================================================================//
template <typename DataType, typename... Args>
DiscreteVariable<DataType> *MultiLevelMeshField::registerCellVariable(
    const std::string &variable_name, Args &&...args)
{
    return registerVariable<DiscreteVariable, DataType>(
        all_cell_variables_, cell_variable_ptrs_, variable_name, total_number_of_cells_,
        std::forward<Args>(args)...);
}
//=============================================================================================//
template <class ExecutionPolicy>
void MultiLevelMeshField::syncCellVariablesToWrite(ExecutionPolicy &ex_policy)
{
    sync_cell_variable_data_(cell_variables_to_write_, ex_policy);
}
//=================================================================================================//
} // namespace SPH
#endif // BASE_MESH_HPP