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
template <typename DataType>
DataType DataValueFromGlobalIndex(MeshVariableData<DataType> *mesh_variable_data,
                                  const Arrayi &global_grid_index,
                                  MeshWithGridDataPackagesType *data_mesh,
                                  size_t *cell_package_index)
{
    constexpr int pkg_size = MeshWithGridDataPackagesType::pkg_size;
    Arrayi cell_index_on_mesh_ = Arrayi::Zero();
    Arrayi local_data_index = Arrayi::Zero();
    for (int n = 0; n != 3; n++)
    {
        size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
        cell_index_on_mesh_[n] = cell_index_in_this_direction;
        local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
    }
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    auto &data = mesh_variable_data[package_index];
    return data[local_data_index[0]][local_data_index[1]][local_data_index[2]];
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP