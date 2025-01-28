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
    auto &data = mesh_variable.Data()[package_index];
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
template <class DataType>
DataType MeshWithGridDataPackages<PKG_SIZE>::
    probeMesh(MeshVariable<DataType> &mesh_variable, const Vecd &position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    size_t package_index = PackageIndexFromCellIndex(cell_index);
    return isInnerDataPackage(cell_index) ? probeDataPackage<DataType>(mesh_variable, package_index, cell_index, position)
                                          : mesh_variable.Data()[package_index][0][0][0];
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
    auto &pkg_data = mesh_variable.Data()[package_index];
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
template <typename InDataType, typename OutDataType>
void MeshWithGridDataPackages<PKG_SIZE>::
    computeGradient(MeshVariable<InDataType> &in_variable,
                    MeshVariable<OutDataType> &out_variable,
                    const size_t package_index)
{
    auto in_variable_data = in_variable.Data();
    auto out_variable_data = out_variable.Data();

    auto &neighborhood = cell_neighborhood_[package_index];
    auto &pkg_data = out_variable_data[package_index];

    for_each_cell_data(
        [&](int i, int j, int k)
        {
            std::pair<size_t, Arrayi> x1 = NeighbourIndexShift(Arrayi(i + 1, j, k), neighborhood);
            std::pair<size_t, Arrayi> x2 = NeighbourIndexShift(Arrayi(i - 1, j, k), neighborhood);
            std::pair<size_t, Arrayi> y1 = NeighbourIndexShift(Arrayi(i, j + 1, k), neighborhood);
            std::pair<size_t, Arrayi> y2 = NeighbourIndexShift(Arrayi(i, j - 1, k), neighborhood);
            std::pair<size_t, Arrayi> z1 = NeighbourIndexShift(Arrayi(i, j, k + 1), neighborhood);
            std::pair<size_t, Arrayi> z2 = NeighbourIndexShift(Arrayi(i, j, k - 1), neighborhood);
            Real dphidx = (in_variable_data[x1.first][x1.second[0]][x1.second[1]][x1.second[2]] -
                           in_variable_data[x2.first][x2.second[0]][x2.second[1]][x2.second[2]]);
            Real dphidy = (in_variable_data[y1.first][y1.second[0]][y1.second[1]][y1.second[2]] -
                           in_variable_data[y2.first][y2.second[0]][y2.second[1]][y2.second[2]]);
            Real dphidz = (in_variable_data[z1.first][z1.second[0]][z1.second[1]][z1.second[2]] -
                           in_variable_data[z2.first][z2.second[0]][z2.second[1]][z2.second[2]]);
            pkg_data[i][j][k] = 0.5 * Vecd(dphidx, dphidy, dphidz) / data_spacing_;
        });
}
//=================================================================================================//
template <int PKG_SIZE>
template <typename DataType>
DataType MeshWithGridDataPackages<PKG_SIZE>::
    CornerAverage(MeshVariable<DataType> &mesh_variable, Arrayi addrs_index, Arrayi corner_direction, CellNeighborhood &neighborhood)
{
    DataType average = ZeroData<DataType>::value;
    auto mesh_variable_data = mesh_variable.Data();
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
template <int PKG_SIZE>
template <class DataType>
DataType MeshWithGridDataPackages<PKG_SIZE>::
    probeDataPackage(MeshVariable<DataType> &mesh_variable, size_t package_index, const Arrayi &cell_index, const Vecd &position)
{
    Arrayi data_index = DataIndexFromPosition(cell_index, position);
    Vecd data_position = DataPositionFromIndex(cell_index, data_index);
    Vecd alpha = (position - data_position) / data_spacing_;
    Vecd beta = Vecd::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    auto mesh_variable_data = mesh_variable.Data();
    NeighbourIndex neighbour_index_1 = NeighbourIndexShift(Arrayi(data_index[0], data_index[1], data_index[2]), neighborhood);
    NeighbourIndex neighbour_index_2 = NeighbourIndexShift(Arrayi(data_index[0] + 1, data_index[1], data_index[2]), neighborhood);
    NeighbourIndex neighbour_index_3 = NeighbourIndexShift(Arrayi(data_index[0], data_index[1] + 1, data_index[2]), neighborhood);
    NeighbourIndex neighbour_index_4 = NeighbourIndexShift(Arrayi(data_index[0] + 1, data_index[1] + 1, data_index[2]), neighborhood);

    DataType bilinear_1 = mesh_variable_data[neighbour_index_1.first][neighbour_index_1.second[0]][neighbour_index_1.second[1]][neighbour_index_1.second[2]] * beta[0] * beta[1] +
                          mesh_variable_data[neighbour_index_2.first][neighbour_index_2.second[0]][neighbour_index_2.second[1]][neighbour_index_2.second[2]] * alpha[0] * beta[1] +
                          mesh_variable_data[neighbour_index_3.first][neighbour_index_3.second[0]][neighbour_index_3.second[1]][neighbour_index_3.second[2]] * beta[0] * alpha[1] +
                          mesh_variable_data[neighbour_index_4.first][neighbour_index_4.second[0]][neighbour_index_4.second[1]][neighbour_index_4.second[2]] * alpha[0] * alpha[1];

    neighbour_index_1 = NeighbourIndexShift(Arrayi(data_index[0], data_index[1], data_index[2] + 1), neighborhood);
    neighbour_index_2 = NeighbourIndexShift(Arrayi(data_index[0] + 1, data_index[1], data_index[2] + 1), neighborhood);
    neighbour_index_3 = NeighbourIndexShift(Arrayi(data_index[0], data_index[1] + 1, data_index[2] + 1), neighborhood);
    neighbour_index_4 = NeighbourIndexShift(Arrayi(data_index[0] + 1, data_index[1] + 1, data_index[2] + 1), neighborhood);

    DataType bilinear_2 = mesh_variable_data[neighbour_index_1.first][neighbour_index_1.second[0]][neighbour_index_1.second[1]][neighbour_index_1.second[2]] * beta[0] * beta[1] +
                          mesh_variable_data[neighbour_index_2.first][neighbour_index_2.second[0]][neighbour_index_2.second[1]][neighbour_index_2.second[2]] * alpha[0] * beta[1] +
                          mesh_variable_data[neighbour_index_3.first][neighbour_index_3.second[0]][neighbour_index_3.second[1]][neighbour_index_3.second[2]] * beta[0] * alpha[1] +
                          mesh_variable_data[neighbour_index_4.first][neighbour_index_4.second[0]][neighbour_index_4.second[1]][neighbour_index_4.second[2]] * alpha[0] * alpha[1];
    return bilinear_1 * beta[2] + bilinear_2 * alpha[2];
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
#endif // MESH_WITH_DATA_PACKAGES_3D_HPP