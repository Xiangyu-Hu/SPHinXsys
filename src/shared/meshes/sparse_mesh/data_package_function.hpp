#ifndef DATA_PACKAGE_FUNCTIONS_HPP
#define DATA_PACKAGE_FUNCTIONS_HPP

#include "data_package_function.h"

namespace SPH
{
//=============================================================================================//
template <int PKG_SIZE>
PackageDataPair NeighbourIndexShift(const Arrayi &shift_index, const CellNeighborhood &neighbour)
{
    PackageDataPair result;
    Arrayi neighbour_index = (shift_index + PKG_SIZE * Arrayi::Ones()) / PKG_SIZE;
    result.first = neighbour(neighbour_index);
    result.second = (shift_index + PKG_SIZE * Arrayi::Ones()) - neighbour_index * PKG_SIZE;
    return result;
}
//=============================================================================================//
template <int PKG_SIZE>
PackageDataPair GeneralNeighbourIndexShift(
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
            PackageDataPair neighbour_index =
                NeighbourIndexShift<PKG_SIZE>(addrs_index + index * corner_direction, neighborhood);
            average += pkg_data[neighbour_index.first](neighbour_index.second);
            count += 1.0;
        });
    return average / count;
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
    PackageDataPair center = NeighbourIndexShift<PKG_SIZE>(data_index, neighborhood);
    PackageDataPair x1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(1, 0), neighborhood);
    PackageDataPair x2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(-1, 0), neighborhood);
    PackageDataPair y1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, 1), neighborhood);
    PackageDataPair y2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array2i(0, -1), neighborhood);
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
    PackageDataPair center = NeighbourIndexShift<PKG_SIZE>(data_index, neighborhood);
    PackageDataPair x1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(1, 0, 0), neighborhood);
    PackageDataPair x2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(-1, 0, 0), neighborhood);
    PackageDataPair y1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 1, 0), neighborhood);
    PackageDataPair y2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, -1, 0), neighborhood);
    PackageDataPair z1 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, 1), neighborhood);
    PackageDataPair z2 = NeighbourIndexShift<PKG_SIZE>(data_index + Array3i(0, 0, -1), neighborhood);
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
} // namespace SPH
#endif // DATA_PACKAGE_FUNCTIONS_HPP