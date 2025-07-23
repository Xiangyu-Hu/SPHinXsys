#ifndef GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP

#include "grid_data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
inline NeighbourIndex NeighbourIndexShift(const Array3i shift_index, const CellNeighborhood3d &neighbour)
{
    NeighbourIndex result;
    Array3i neighbour_index = (shift_index + PKG_SIZE * Array3i::Ones()) / PKG_SIZE;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]][neighbour_index[2]];
    result.second = (shift_index + PKG_SIZE * Array3i::Ones()) - neighbour_index * PKG_SIZE;

    return result;
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType CornerAverage(PackageDataMatrix3d<DataType, PKG_SIZE> *pkg_data, Array3i addrs_index,
                       Array3i corner_direction, const CellNeighborhood3d &neighborhood, DataType zero)
{
    DataType average = zero;
    for (int i = 0; i != 2; ++i)
        for (int j = 0; j != 2; ++j)
            for (int k = 0; k != 2; ++k)
            {
                int x_index = addrs_index[0] + i * corner_direction[0];
                int y_index = addrs_index[1] + j * corner_direction[1];
                int z_index = addrs_index[2] + k * corner_direction[2];
                NeighbourIndex neighbour_index =
                    NeighbourIndexShift<PKG_SIZE>(Array3i(x_index, y_index, z_index), neighborhood);
                average += pkg_data[neighbour_index.first]
                                   [neighbour_index.second[0]]
                                   [neighbour_index.second[1]]
                                   [neighbour_index.second[2]];
            }
    return average * 0.125;
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType DataValueFromGlobalIndex(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data,
                                  const Array3i &global_grid_index,
                                  MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
                                  size_t *cell_package_index)
{
    Array3i cell_index_on_mesh_ = global_grid_index / PKG_SIZE;
    Array3i local_index = global_grid_index - cell_index_on_mesh_ * PKG_SIZE;
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    return pkg_data[package_index][local_index[0]][local_index[1]][local_index[2]];
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP