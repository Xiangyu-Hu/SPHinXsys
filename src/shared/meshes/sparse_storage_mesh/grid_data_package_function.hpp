#ifndef GRID_DATA_PACKAGE_FUNCTIONS_HPP
#define GRID_DATA_PACKAGE_FUNCTIONS_HPP

#include "grid_data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
PackageGridPair NeighbourIndexShift(const Arrayi &shift_index, const CellNeighborhood &neighbour)
{
    PackageGridPair result;
    Arrayi neighbour_index = (shift_index + PKG_SIZE * Arrayi::Ones()) / PKG_SIZE;
    result.first = neighbour(neighbour_index);
    result.second = (shift_index + PKG_SIZE * Arrayi::Ones()) - neighbour_index * PKG_SIZE;
    return result;
}
//=============================================================================================//
template <int PKG_SIZE>
PackageGridPair GeneralNeighbourIndexShift(
    UnsignedInt package_index, CellNeighborhood *neighbour, const Arrayi &shift_index)
{
    Arrayi cell_shift = shift_index / PKG_SIZE;
    Arrayi residual = shift_index - cell_shift * PKG_SIZE;
    while (!cell_shift.isZero())
    {
        for (UnsignedInt i = 0; i != Dimensions; ++i)
            if (cell_shift[i] != 0)
            {
                Arrayi step = Arrayi::Zero();
                step[i] = cell_shift[i] > 0 ? -1 : 1;
                Arrayi neighbour_index = Arrayi::Ones() - step;
                package_index = neighbour[package_index](neighbour_index);
                cell_shift[i] += step[i];
            }
    }
    return NeighbourIndexShift<PKG_SIZE>(residual, neighbour[package_index]);
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType CornerAverage(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data, Arrayi addrs_index,
                       Arrayi corner_direction, const CellNeighborhood &neighborhood, DataType zero)
{
    DataType average = zero;
    Real count = 0.0;
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Ones() * 2,
        [&](const Arrayi &index)
        {
            PackageGridPair neighbour_index =
                NeighbourIndexShift<PKG_SIZE>(addrs_index + index * corner_direction, neighborhood);
            average += pkg_data[neighbour_index.first](neighbour_index.second);
            count += 1.0;
        });
    return average / count;
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType DataValueFromGlobalIndex(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data,
                                  const Arrayi &global_grid_index,
                                  MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
                                  size_t *cell_package_index)
{
    Arrayi cell_index_on_mesh_ = global_grid_index / PKG_SIZE;
    Arrayi local_index = global_grid_index - cell_index_on_mesh_ * PKG_SIZE;
    size_t package_index = data_mesh->PackageIndexFromCellIndex(cell_package_index, cell_index_on_mesh_);
    return pkg_data[package_index](local_index);
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
template <class ExecutionPolicy>
ProbeMesh<DataType, PKG_SIZE>::ProbeMesh(
    const ExecutionPolicy &ex_policy, MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
    const std::string variable_name)
    : pkg_data_(data_mesh->template getMeshVariable<DataType>(variable_name)->DelegatedData(ex_policy)),
      index_handler_(data_mesh->index_handler_.DelegatedData(ex_policy)),
      cell_pkg_index_(data_mesh->cell_pkg_index_.DelegatedData(ex_policy)),
      cell_neighborhood_(data_mesh->cell_neighborhood_.DelegatedData(ex_policy)) {}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::operator()(const Vecd &position)
{
    Arrayi cell_index = index_handler_->CellIndexFromPosition(position);
    size_t package_index = index_handler_->PackageIndexFromCellIndex(cell_pkg_index_, cell_index);
    return package_index > 1 ? probeDataPackage(package_index, cell_index, position)
                             : pkg_data_[package_index](Arrayi::Zero());
}
//=============================================================================================//
template <typename CellDataType, typename PackageDataType, size_t PKG_SIZE, typename FunctionByGrid>
CellDataType assignByGrid(PackageDataMatrix<PackageDataType, PKG_SIZE> &pkg_data,
                          const FunctionByGrid &function_by_grid, CellDataType inital_value)
{
    CellDataType value = inital_value;
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Ones() * PKG_SIZE,
        [&](const Arrayi &index)
        {
            value = function_by_grid(pkg_data(index), index, value);
        });
    return value;
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::probeDataPackage(
    size_t package_index, const Array2i &cell_index, const Vec2d &position)
{
    Array2i data_index = index_handler_->DataIndexFromPosition(cell_index, position);
    Vec2d data_position = index_handler_->DataPositionFromIndex(cell_index, data_index);
    Vec2d alpha = (position - data_position) / index_handler_->data_spacing_;
    Vec2d beta = Vec2d::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    PackageGridPair neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, 0), neighborhood);
    PackageGridPair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(1, 0), neighborhood);
    PackageGridPair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, 1), neighborhood);
    PackageGridPair neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(1, 1), neighborhood);

    return pkg_data_[neighbour_index_1.first](neighbour_index_1.second) * beta[0] * beta[1] +
           pkg_data_[neighbour_index_2.first](neighbour_index_2.second) * alpha[0] * beta[1] +
           pkg_data_[neighbour_index_3.first](neighbour_index_3.second) * beta[0] * alpha[1] +
           pkg_data_[neighbour_index_4.first](neighbour_index_4.second) * alpha[0] * alpha[1];
}
//=============================================================================================//
template <typename DataType, size_t PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::probeDataPackage(
    size_t package_index, const Array3i &cell_index, const Vec3d &position)
{
    Array3i data_index = index_handler_->DataIndexFromPosition(cell_index, position);
    Vec3d data_position = index_handler_->DataPositionFromIndex(cell_index, data_index);
    Vec3d alpha = (position - data_position) / index_handler_->data_spacing_;
    Vec3d beta = Vec3d::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    PackageGridPair neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 0), neighborhood);
    PackageGridPair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 0), neighborhood);
    PackageGridPair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 0), neighborhood);
    PackageGridPair neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 1, 0), neighborhood);

    DataType bilinear_1 =
        pkg_data_[neighbour_index_1.first](neighbour_index_1.second) * beta[0] * beta[1] +
        pkg_data_[neighbour_index_2.first](neighbour_index_2.second) * alpha[0] * beta[1] +
        pkg_data_[neighbour_index_3.first](neighbour_index_3.second) * beta[0] * alpha[1] +
        pkg_data_[neighbour_index_4.first](neighbour_index_4.second) * alpha[0] * alpha[1];

    neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 1), neighborhood);
    neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 1), neighborhood);
    neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 1), neighborhood);
    neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 1, 1), neighborhood);

    DataType bilinear_2 =
        pkg_data_[neighbour_index_1.first](neighbour_index_1.second) * beta[0] * beta[1] +
        pkg_data_[neighbour_index_2.first](neighbour_index_2.second) * alpha[0] * beta[1] +
        pkg_data_[neighbour_index_3.first](neighbour_index_3.second) * beta[0] * alpha[1] +
        pkg_data_[neighbour_index_4.first](neighbour_index_4.second) * alpha[0] * alpha[1];

    return bilinear_1 * beta[2] + bilinear_2 * alpha[2];
}
//=============================================================================================//
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_HPP