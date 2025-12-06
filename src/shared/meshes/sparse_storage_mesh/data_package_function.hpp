#ifndef DATA_PACKAGE_FUNCTIONS_HPP
#define DATA_PACKAGE_FUNCTIONS_HPP

#include "data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
DataPackagePair NeighbourIndexShift(const Arrayi &shift_index, const CellNeighborhood &neighbour)
{
    DataPackagePair result;
    Arrayi neighbour_index = (shift_index + PKG_SIZE * Arrayi::Ones()) / PKG_SIZE;
    result.first = neighbour(neighbour_index);
    result.second = (shift_index + PKG_SIZE * Arrayi::Ones()) - neighbour_index * PKG_SIZE;
    return result;
}
//=============================================================================================//
template <int PKG_SIZE>
DataPackagePair GeneralNeighbourIndexShift(
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
template <typename DataType, int PKG_SIZE>
DataType CornerAverage(PackageData<DataType, PKG_SIZE> *pkg_data, Arrayi addrs_index,
                       Arrayi corner_direction, const CellNeighborhood &neighborhood, DataType zero)
{
    DataType average = zero;
    Real count = 0.0;
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Ones() * 2,
        [&](const Arrayi &index)
        {
            DataPackagePair neighbour_index =
                NeighbourIndexShift<PKG_SIZE>(addrs_index + index * corner_direction, neighborhood);
            average += pkg_data[neighbour_index.first](neighbour_index.second);
            count += 1.0;
        });
    return average / count;
}
//=============================================================================================//
template <typename DataType, int PKG_SIZE>
template <class ExecutionPolicy>
ProbeMesh<DataType, PKG_SIZE>::ProbeMesh(
    const ExecutionPolicy &ex_policy, MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
    const std::string variable_name)
    : pkg_data_(data_mesh->template getMeshVariable<DataType>(variable_name)->DelegatedData(ex_policy)),
      index_handler_(data_mesh->getIndexHandler()),
      cell_pkg_index_(data_mesh->getCellPackageIndex().DelegatedData(ex_policy)),
      cell_neighborhood_(data_mesh->getCellNeighborhood().DelegatedData(ex_policy)) {}
//=============================================================================================//
template <typename DataType, int PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::operator()(const Vecd &position)
{
    Arrayi cell_index = index_handler_.CellIndexFromPosition(position);
    UnsignedInt package_index = index_handler_.PackageIndexFromCellIndex(cell_pkg_index_, cell_index);
    return package_index > 1 ? probeDataPackage(package_index, cell_index, position)
                             : pkg_data_[package_index](Arrayi::Zero());
}
//=============================================================================================//
template <typename DataType, int PKG_SIZE, typename FunctionByIndex>
void assignByDataIndex(PackageData<DataType, PKG_SIZE> &pkg_data, const FunctionByIndex &function_by_index)
{
    mesh_for_each(
        Arrayi::Zero(), Arrayi::Ones() * PKG_SIZE,
        [&](const Arrayi &data_index)
        {
            pkg_data(data_index) = function_by_index(data_index);
        });
}
//=============================================================================================//
template <int PKG_SIZE, typename RegularizeFunction>
Vec2d regularizedCentralDifference(
    PackageData<Real, PKG_SIZE> *input, const CellNeighborhood2d &neighborhood,
    const Array2i &data_index, const RegularizeFunction &regularize_function)
{
    DataPackagePair center = NeighbourIndexShift<PKG_SIZE>(data_index, neighborhood);
    DataPackagePair x1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(1, 0), neighborhood);
    DataPackagePair x2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(-1, 0), neighborhood);
    DataPackagePair y1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, 1), neighborhood);
    DataPackagePair y2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, -1), neighborhood);
    Real dphidx_p = input[x1.first](x1.second) - input[center.first](center.second);
    Real dphidx_m = input[center.first](center.second) - input[x2.first](x2.second);
    Real dphidx = regularize_function(dphidx_p, dphidx_m);
    Real dphidy_p = input[y1.first](y1.second) - input[center.first](center.second);
    Real dphidy_m = input[center.first](center.second) - input[y2.first](y2.second);
    Real dphidy = regularize_function(dphidy_p, dphidy_m);
    return Vec2d(dphidx, dphidy);
}
//=============================================================================================//
template <int PKG_SIZE, typename RegularizeFunction>
Vec3d regularizedCentralDifference(
    PackageData<Real, PKG_SIZE> *input, const CellNeighborhood3d &neighborhood,
    const Array3i &data_index, const RegularizeFunction &regularize_function)
{
    DataPackagePair center = NeighbourIndexShift<PKG_SIZE>(data_index, neighborhood);
    DataPackagePair x1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 0), neighborhood);
    DataPackagePair x2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(-1, 0, 0), neighborhood);
    DataPackagePair y1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 0), neighborhood);
    DataPackagePair y2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, -1, 0), neighborhood);
    DataPackagePair z1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 1), neighborhood);
    DataPackagePair z2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, -1), neighborhood);
    Real dphidx_p = input[x1.first](x1.second) - input[center.first](center.second);
    Real dphidx_m = input[center.first](center.second) - input[x2.first](x2.second);
    Real dphidx = regularize_function(dphidx_p, dphidx_m);
    Real dphidy_p = input[y1.first](y1.second) - input[center.first](center.second);
    Real dphidy_m = input[center.first](center.second) - input[y2.first](y2.second);
    Real dphidy = regularize_function(dphidy_p, dphidy_m);
    Real dphidz_p = input[z1.first](z1.second) - input[center.first](center.second);
    Real dphidz_m = input[center.first](center.second) - input[z2.first](z2.second);
    Real dphidz = regularize_function(dphidz_p, dphidz_m);
    return Vec3d(dphidx, dphidy, dphidz);
}
//=============================================================================================//
template <typename DataType, int PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::probeDataPackage(
    UnsignedInt package_index, const Array2i &cell_index, const Vec2d &position)
{
    Array2i data_index = index_handler_.DataIndexFromPosition(cell_index, position);
    Vec2d data_position = index_handler_.DataPositionFromIndex(cell_index, data_index);
    Vec2d alpha = (position - data_position) / index_handler_.DataSpacing();
    Vec2d beta = Vec2d::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    DataPackagePair neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, 0), neighborhood);
    DataPackagePair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(1, 0), neighborhood);
    DataPackagePair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, 1), neighborhood);
    DataPackagePair neighbour_index_4 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(1, 1), neighborhood);

    return pkg_data_[neighbour_index_1.first](neighbour_index_1.second) * beta[0] * beta[1] +
           pkg_data_[neighbour_index_2.first](neighbour_index_2.second) * alpha[0] * beta[1] +
           pkg_data_[neighbour_index_3.first](neighbour_index_3.second) * beta[0] * alpha[1] +
           pkg_data_[neighbour_index_4.first](neighbour_index_4.second) * alpha[0] * alpha[1];
}
//=============================================================================================//
template <typename DataType, int PKG_SIZE>
DataType ProbeMesh<DataType, PKG_SIZE>::probeDataPackage(
    UnsignedInt package_index, const Array3i &cell_index, const Vec3d &position)
{
    Array3i data_index = index_handler_.DataIndexFromPosition(cell_index, position);
    Vec3d data_position = index_handler_.DataPositionFromIndex(cell_index, data_index);
    Vec3d alpha = (position - data_position) / index_handler_.DataSpacing();
    Vec3d beta = Vec3d::Ones() - alpha;

    auto &neighborhood = cell_neighborhood_[package_index];
    DataPackagePair neighbour_index_1 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 0), neighborhood);
    DataPackagePair neighbour_index_2 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 0), neighborhood);
    DataPackagePair neighbour_index_3 =
        NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 0), neighborhood);
    DataPackagePair neighbour_index_4 =
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
#endif // DATA_PACKAGE_FUNCTIONS_HPP