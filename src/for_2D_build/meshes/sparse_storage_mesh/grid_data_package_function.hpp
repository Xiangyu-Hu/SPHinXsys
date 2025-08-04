#ifndef GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP

#include "grid_data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
PackageGridPair NeighbourIndexShift(const Array2i &shift_index, const CellNeighborhood2d &neighbour)
{
    PackageGridPair result;
    Array2i neighbour_index = (shift_index + PKG_SIZE * Array2i::Ones()) / PKG_SIZE;
    result.first = neighbour[neighbour_index[0]][neighbour_index[1]];
    result.second = (shift_index + PKG_SIZE * Array2i::Ones()) - neighbour_index * PKG_SIZE;

    return result;
}
//=============================================================================================//
template <int PKG_SIZE>
PackageGridPair GeneralNeighbourIndexShift(
    UnsignedInt package_index, CellNeighborhood *neighbour, const Array2i &shift_index)
{
    Array2i cell_shift = shift_index / PKG_SIZE;
    for (UnsignedInt i = 0; i != Dimensions; ++i)
    {
        int n = cell_shift[i];
        Array2i step = Array2i::Zero();
        step[i] = n > 0 ? -1 : 1;
        Array2i neighbour_index = Array2i::Ones() - step;
        for (int j = n; j != 0; j += step[i])
        {
            package_index = neighbour[package_index]
                                     [neighbour_index[0]]
                                     [neighbour_index[1]];
        }
    }
    Array2i residual = shift_index - cell_shift * PKG_SIZE;
    return NeighbourIndexShift<PKG_SIZE>(residual, neighbour[package_index]);
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
            PackageGridPair neighbour_index =
                NeighbourIndexShift<PKG_SIZE>(Array2i(x_index, y_index), neighborhood);
            average += pkg_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]];
        }
    return average * 0.25;
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType DataValueFromGlobalIndex(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data,
                                  const Array2i &global_grid_index,
                                  MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
                                  size_t *cell_package_index)
{
    Array2i cell_index_on_mesh_ = global_grid_index / PKG_SIZE;
    Array2i local_index = global_grid_index - cell_index_on_mesh_ * PKG_SIZE;
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    return pkg_data[package_index][local_index[0]][local_index[1]];
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::operator()(const Vecd &position)
{
    Arrayi cell_index = index_handler_->CellIndexFromPosition(position);
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_package_index_, cell_index);
    return package_index > 1 ? probeDataPackage(package_index, cell_index, position)
                             : pkg_data_[package_index][0][0];
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
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(0, 0), neighborhood);
    PackageGridPair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(1, 0), neighborhood);
    PackageGridPair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(0, 1), neighborhood);
    PackageGridPair neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Arrayi(1, 1), neighborhood);

    return pkg_data_[neighbour_index_1.first][neighbour_index_1.second[0]][neighbour_index_1.second[1]] * beta[0] * beta[1] +
           pkg_data_[neighbour_index_2.first][neighbour_index_2.second[0]][neighbour_index_2.second[1]] * alpha[0] * beta[1] +
           pkg_data_[neighbour_index_3.first][neighbour_index_3.second[0]][neighbour_index_3.second[1]] * beta[0] * alpha[1] +
           pkg_data_[neighbour_index_4.first][neighbour_index_4.second[0]][neighbour_index_4.second[1]] * alpha[0] * alpha[1];
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_2D_HPP