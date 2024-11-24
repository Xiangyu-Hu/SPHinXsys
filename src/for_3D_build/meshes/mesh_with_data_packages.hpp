/**
 * @file    mesh_with_data_packages.hpp
 * @brief   Implementation for 3d builds.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_3D_HPP
#define MESH_WITH_DATA_PACKAGES_3D_HPP

#include "mesh_with_data_packages.h"
#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType MeshWithGridDataPackages<PKG_SIZE>::
    DataValueFromGlobalIndex(MeshVariable<DataType> &mesh_variable,
                             const Arrayi &global_grid_index)
{
    Arrayi cell_index_on_mesh_ = Arrayi::Zero();
    Arrayi local_data_index = Arrayi::Zero();
    for (int n = 0; n != 3; n++)
    {
        size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
        cell_index_on_mesh_[n] = cell_index_in_this_direction;
        local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
    }
    size_t package_index = PackageIndexFromCellIndex(cell_index_on_mesh_);
    auto &data = mesh_variable.DataField()[package_index];
    return data[local_data_index[0]][local_data_index[1]][local_data_index[2]];
}
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::allocateIndexDataMatrix()
{
    Allocate3dArray(index_data_mesh_, all_cells_);
}
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::deleteIndexDataMatrix()
{
    Delete3dArray(index_data_mesh_, all_cells_);
}
//=================================================================================================//
template <int PKG_SIZE>
template <typename FunctionOnData>
void MeshWithGridDataPackages<PKG_SIZE>::
    for_each_cell_data(const FunctionOnData &function)
{
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
            for (int k = 0; k != pkg_size; ++k)
            {
                function(i, j, k);
            }
}
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::
    assignDataPackageIndex(const Arrayi &cell_index, const size_t package_index)
{
    index_data_mesh_[cell_index[0]][cell_index[1]][cell_index[2]] = package_index;
}
//=================================================================================================//
template <int PKG_SIZE>
size_t MeshWithGridDataPackages<PKG_SIZE>::
    PackageIndexFromCellIndex(const Arrayi &cell_index)
{
    return index_data_mesh_[cell_index[0]][cell_index[1]][cell_index[2]];
}
//=================================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::
    isSingularDataPackage(const Arrayi &cell_index)
{
    return index_data_mesh_[cell_index[0]][cell_index[1]][cell_index[2]] < 2;
}
//=================================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::
    isInnerDataPackage(const Arrayi &cell_index)
{
    return index_data_mesh_[cell_index[0]][cell_index[1]][cell_index[2]] > 1;
}
//=================================================================================================//
template <int PKG_SIZE>
std::pair<size_t, Arrayi> MeshWithGridDataPackages<PKG_SIZE>::
    NeighbourIndexShift(const Arrayi shift_index, const CellNeighborhood &neighbour)
{
    std::pair<size_t, Arrayi> result;
    Arrayi neighbour_index = (shift_index + pkg_size * Arrayi::Ones()) / pkg_size;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]][neighbour_index[2]];
    result.second = (shift_index + pkg_size * Arrayi::Ones()) - neighbour_index * pkg_size;

    return result;
}
//=================================================================================================//
template <int PKG_SIZE>
template <typename DataType, typename FunctionByPosition>
void MeshWithGridDataPackages<PKG_SIZE>::
    assignByPosition(MeshVariable<DataType> &mesh_variable,
                     const Arrayi &cell_index,
                     const FunctionByPosition &function_by_position)
{
    size_t package_index = PackageIndexFromCellIndex(cell_index);
    auto &pkg_data = mesh_variable.DataField()[package_index];
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
            for (int k = 0; k != pkg_size; ++k)
            {
                Vec3d position = DataPositionFromIndex(cell_index, Arrayi(i, j, k));
                pkg_data[i][j][k] = function_by_position(position);
            }
}
//=================================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType MeshWithGridDataPackages<PKG_SIZE>::
    CornerAverage(MeshVariable<DataType> &mesh_variable, Arrayi addrs_index, Arrayi corner_direction, CellNeighborhood &neighborhood)
{
    DataType average = ZeroData<DataType>::value;
    auto mesh_variable_data = mesh_variable.DataField();
    for (int i = 0; i != 2; ++i)
        for (int j = 0; j != 2; ++j)
            for (int k = 0; k != 2; ++k)
            {
                int x_index = addrs_index[0] + i * corner_direction[0];
                int y_index = addrs_index[1] + j * corner_direction[1];
                int z_index = addrs_index[2] + k * corner_direction[2];
                NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(x_index, y_index, z_index), neighborhood);
                average += mesh_variable_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
            }
    return average * 0.125;
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
#endif // MESH_WITH_DATA_PACKAGES_3D_HPP