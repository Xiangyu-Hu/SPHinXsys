#ifndef GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP

#include "grid_data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
NeighbourIndex NeighbourIndexShift(const Array2i shift_index, const CellNeighborhood2d &neighbour)
{
    NeighbourIndex result;
    Array2i neighbour_index = (shift_index + PKG_SIZE * Array2i::Ones()) / PKG_SIZE;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]];
    result.second = (shift_index + PKG_SIZE * Array2i::Ones()) - neighbour_index * PKG_SIZE;

    return result;
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType CornerAverage(PackageDataMatrix2d<DataType, PKG_SIZE> *pkg_data, Array2i addrs_index,
                       Array2i corner_direction, const CellNeighborhood2d &neighborhood, DataType zero)
{
    DataType average = ZeroData<DataType>::value;
    for (int i = 0; i != 2; ++i)
        for (int j = 0; j != 2; ++j)
        {
            int x_index = addrs_index[0] + i * corner_direction[0];
            int y_index = addrs_index[1] + j * corner_direction[1];
            NeighbourIndex neighbour_index =
                NeighbourIndexShift<PKG_SIZE>(Array2i(x_index, y_index), neighborhood);
            average += pkg_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
        }
    return average * 0.25;
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
    for (int n = 0; n != 2; n++)
    {
        size_t cell_index_in_this_direction = global_grid_index[n] / pkg_size;
        cell_index_on_mesh_[n] = cell_index_in_this_direction;
        local_data_index[n] = global_grid_index[n] - cell_index_in_this_direction * pkg_size;
    }
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    auto &data = mesh_variable_data[package_index];
    return data[local_data_index[0]][local_data_index[1]];
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP