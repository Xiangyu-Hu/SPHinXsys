#ifndef GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP

#include "grid_data_package_function.h"

namespace SPH
{
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
                PackageGridPair neighbour_index =
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
template <typename DataType, size_t PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::operator()(const Vecd &position)
{
    Arrayi cell_index = index_handler_->CellIndexFromPosition(position);
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_package_index_, cell_index);
    return package_index > 1 ? probeDataPackage(package_index, cell_index, position)
                             : pkg_data_[package_index][0][0][0];
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::probeDataPackage(
    size_t package_index, const Arrayi &cell_index, const Vecd &position)
{
    Arrayi data_index = index_handler_->DataIndexFromPosition(cell_index, position);
    Vecd data_position = index_handler_->DataPositionFromIndex(cell_index, data_index);
    Vecd alpha = (position - data_position) / index_handler_->data_spacing_;
    Vecd beta = Vecd::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    PackageGridPair neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(0, 0, 0), neighborhood);
    PackageGridPair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(1, 0, 0), neighborhood);
    PackageGridPair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(0, 1, 0), neighborhood);
    PackageGridPair neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(1, 1, 0), neighborhood);

    DataType bilinear_1 =
        pkg_data_[neighbour_index_1.first]
                 [neighbour_index_1.second[0]]
                 [neighbour_index_1.second[1]]
                 [neighbour_index_1.second[2]] *
            beta[0] * beta[1] +
        pkg_data_[neighbour_index_2.first]
                 [neighbour_index_2.second[0]]
                 [neighbour_index_2.second[1]]
                 [neighbour_index_2.second[2]] *
            alpha[0] * beta[1] +
        pkg_data_[neighbour_index_3.first]
                 [neighbour_index_3.second[0]]
                 [neighbour_index_3.second[1]]
                 [neighbour_index_3.second[2]] *
            beta[0] * alpha[1] +
        pkg_data_[neighbour_index_4.first]
                 [neighbour_index_4.second[0]]
                 [neighbour_index_4.second[1]]
                 [neighbour_index_4.second[2]] *
            alpha[0] * alpha[1];

    neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(0, 0, 1), neighborhood);
    neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(1, 0, 1), neighborhood);
    neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(0, 1, 1), neighborhood);
    neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(1, 1, 1), neighborhood);

    DataType bilinear_2 =
        pkg_data_[neighbour_index_1.first]
                 [neighbour_index_1.second[0]]
                 [neighbour_index_1.second[1]]
                 [neighbour_index_1.second[2]] *
            beta[0] * beta[1] +
        pkg_data_[neighbour_index_2.first]
                 [neighbour_index_2.second[0]]
                 [neighbour_index_2.second[1]]
                 [neighbour_index_2.second[2]] *
            alpha[0] * beta[1] +
        pkg_data_[neighbour_index_3.first]
                 [neighbour_index_3.second[0]]
                 [neighbour_index_3.second[1]]
                 [neighbour_index_3.second[2]] *
            beta[0] * alpha[1] +
        pkg_data_[neighbour_index_4.first]
                 [neighbour_index_4.second[0]]
                 [neighbour_index_4.second[1]]
                 [neighbour_index_4.second[2]] *
            alpha[0] * alpha[1];

    return bilinear_1 * beta[2] + bilinear_2 * alpha[2];
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_3D_HPP