/**
 * @file    mesh_with_data_packages.hpp
 * @brief   Implementation for 3d builds.
 * @author  Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_3D_HPP
#define MESH_WITH_DATA_PACKAGES_3D_HPP

#include "mesh_with_data_packages.h"

namespace SPH
{
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename FunctionOnData>
void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
    for_each_data(const FunctionOnData &function)
{
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
            for (int k = 0; k != pkg_size; ++k)
            {
                function(i, j, k);
            }
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename FunctionOnAddress>
void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
    for_each_addrs(const FunctionOnAddress &function)
{
    for (int i = pkg_addrs_buffer; i != pkg_ops_end; ++i)
        for (int j = pkg_addrs_buffer; j != pkg_ops_end; ++j)
            for (int k = pkg_addrs_buffer; k != pkg_ops_end; ++k)
            {
                function(i, j, k);
            }
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <class DataType>
DataType GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
    probeDataPackage(PackageDataAddress<DataType> &pkg_data_addrs, const Vecd &position)
{
    Arrayi grid_idx = CellIndexFromPosition(position);
    Vecd grid_pos = GridPositionFromIndex(grid_idx);
    Vecd alpha = (position - grid_pos) / grid_spacing_;
    Vecd beta = Vecd::Ones() - alpha;

    DataType bilinear_1 = *pkg_data_addrs[grid_idx[0]][grid_idx[1]][grid_idx[2]] * beta[0] * beta[1] +
                          *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]][grid_idx[2]] * alpha[0] * beta[1] +
                          *pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1][grid_idx[2]] * beta[0] * alpha[1] +
                          *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1][grid_idx[2]] * alpha[0] * alpha[1];
    DataType bilinear_2 = *pkg_data_addrs[grid_idx[0]][grid_idx[1]][grid_idx[2] + 1] * beta[0] * beta[1] +
                          *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]][grid_idx[2] + 1] * alpha[0] * beta[1] +
                          *pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1][grid_idx[2] + 1] * beta[0] * alpha[1] +
                          *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1][grid_idx[2] + 1] * alpha[0] * alpha[1];
    return bilinear_1 * beta[2] + bilinear_2 * alpha[2];
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename InDataType, typename OutDataType>
void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
    computeGradient(const MeshVariable<InDataType> &in_variable,
                    const MeshVariable<OutDataType> &out_variable)
{
    auto &in_variable_addrs = getPackageDataAddress(in_variable);
    auto &out_variable_addrs = getPackageDataAddress(out_variable);

    for_each_addrs(
        [&](int i, int j, int k)
        {
            Real dphidx = (*in_variable_addrs[i + 1][j][k] - *in_variable_addrs[i - 1][j][k]);
            Real dphidy = (*in_variable_addrs[i][j + 1][k] - *in_variable_addrs[i][j - 1][k]);
            Real dphidz = (*in_variable_addrs[i][j][k + 1] - *in_variable_addrs[i][j][k - 1]);
            *out_variable_addrs[i][j][k] = 0.5 * Vecd(dphidx, dphidy, dphidz) / grid_spacing_;
        });
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename DataType, typename FunctionByPosition>
void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
    assignByPosition(const MeshVariable<DataType> &mesh_variable,
                     const FunctionByPosition &function_by_position)
{
    auto &pkg_data = getPackageData(mesh_variable);
    for (int i = 0; i != pkg_size; ++i)
        for (int j = 0; j != pkg_size; ++j)
            for (int k = 0; k != pkg_size; ++k)
            {
                Vec3d position = DataPositionFromIndex(Vec3d(i, j, k));
                pkg_data[i][j][k] = function_by_position(position);
            }
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename DataType>
void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::AssignSingularPackageDataAddress<DataType>::
operator()(DataContainerAssemble<PackageData> &all_pkg_data,
           DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    for (size_t l = 0; l != std::get<type_index>(all_pkg_data).size(); ++l)
    {
        PackageData<DataType> &pkg_data = std::get<type_index>(all_pkg_data)[l];
        PackageDataAddress<DataType> &pkg_data_addrs = std::get<type_index>(all_pkg_data_addrs)[l];
        for (int i = 0; i != pkg_addrs_size; ++i)
            for (int j = 0; j != pkg_addrs_size; ++j)
                for (int k = 0; k != pkg_addrs_size; ++k)
                {
                    pkg_data_addrs[i][j][k] = &pkg_data[0][0][0];
                }
    }
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename DataType>
void GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::AssignPackageDataAddress<DataType>::
operator()(DataContainerAssemble<PackageDataAddress> &all_pkg_data_addrs,
           const Arrayi &addrs_index,
           DataContainerAssemble<PackageData> &all_pkg_data,
           const Arrayi &data_index)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    for (size_t l = 0; l != std::get<type_index>(all_pkg_data).size(); ++l)
    {
        PackageData<DataType> &pkg_data = std::get<type_index>(all_pkg_data)[l];
        PackageDataAddress<DataType> &pkg_data_addrs = std::get<type_index>(all_pkg_data_addrs)[l];
        pkg_data_addrs[addrs_index[0]][addrs_index[1]][addrs_index[2]] =
            &pkg_data[data_index[0]][data_index[1]][data_index[2]];
    }
}
//=================================================================================================//
template <int PKG_SIZE, int ADDRS_BUFFER>
template <typename DataType>
DataType GridDataPackage<PKG_SIZE, ADDRS_BUFFER>::
    CornerAverage(PackageDataAddress<DataType> &pkg_data_addrs, Arrayi addrs_index, Arrayi corner_direction)
{
    DataType average = ZeroData<DataType>::value;
    for (int i = 0; i != 2; ++i)
        for (int j = 0; j != 2; ++j)
            for (int k = 0; k != 2; ++k)
            {
                int x_index = addrs_index[0] + i * corner_direction[0];
                int y_index = addrs_index[1] + j * corner_direction[1];
                int z_index = addrs_index[2] + k * corner_direction[2];
                average += *pkg_data_addrs[x_index][y_index][z_index];
            }
    return average * 0.125;
}
//=================================================================================================//
template <class GridDataPackageType>
template <typename DataType>
DataType MeshWithGridDataPackages<GridDataPackageType>::
    DataValueFromGlobalIndex(const MeshVariable<DataType> &mesh_variable,
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
    auto &data = data_pkg_addrs_[cell_index_on_mesh_[0]][cell_index_on_mesh_[1]][cell_index_on_mesh_[2]]
                     ->getPackageData(mesh_variable);
    return data[local_data_index[0]][local_data_index[1]][local_data_index[2]];
}
//=================================================================================================//
template <class GridDataPackageType>
void MeshWithGridDataPackages<GridDataPackageType>::
    initializePackageAddressesInACell(const Arrayi &cell_index)
{
    int i = cell_index[0];
    int j = cell_index[1];
    int k = cell_index[2];

    GridDataPackageType *data_pkg = data_pkg_addrs_[i][j][k];
    if (data_pkg->isInnerPackage())
    {
        for (int l = 0; l != pkg_addrs_size; ++l)
            for (int m = 0; m != pkg_addrs_size; ++m)
                for (int n = 0; n != pkg_addrs_size; ++n)
                {
                    std::pair<int, int> x_pair = CellShiftAndDataIndex(l);
                    std::pair<int, int> y_pair = CellShiftAndDataIndex(m);
                    std::pair<int, int> z_pair = CellShiftAndDataIndex(n);

                    data_pkg->assignPackageDataAddress(
                        Arrayi(l, m, n),
                        data_pkg_addrs_[i + x_pair.first][j + y_pair.first][k + z_pair.first],
                        Arrayi(x_pair.second, y_pair.second, z_pair.second));
                }
    }
}
//=================================================================================================//
template <class GridDataPackageType>
void MeshWithGridDataPackages<GridDataPackageType>::allocateMeshDataMatrix()
{
    Allocate3dArray(data_pkg_addrs_, all_cells_);
}
//=================================================================================================//
template <class GridDataPackageType>
void MeshWithGridDataPackages<GridDataPackageType>::deleteMeshDataMatrix()
{
    Delete3dArray(data_pkg_addrs_, all_cells_);
}
//=================================================================================================//
template <class GridDataPackageType>
void MeshWithGridDataPackages<GridDataPackageType>::
    assignDataPackageAddress(const Arrayi &cell_index, GridDataPackageType *data_pkg)
{
    data_pkg_addrs_[cell_index[0]][cell_index[1]][cell_index[2]] = data_pkg;
}
//=================================================================================================//
template <class GridDataPackageType>
GridDataPackageType *MeshWithGridDataPackages<GridDataPackageType>::
    DataPackageFromCellIndex(const Arrayi &cell_index)
{
    return data_pkg_addrs_[cell_index[0]][cell_index[1]][cell_index[2]];
}
//=================================================================================================//
template <class GridDataPackageType>
template <class DataType>
DataType MeshWithGridDataPackages<GridDataPackageType>::
    probeMesh(const MeshVariable<DataType> &mesh_variable, const Vecd &position)
{
    Arrayi index = CellIndexFromPosition(position);
    GridDataPackageType *data_pkg = data_pkg_addrs_[index[0]][index[1]][index[2]];
    auto &pkg_data_addrs = data_pkg->getPackageDataAddress(mesh_variable);
    return data_pkg->isInnerPackage() ? data_pkg->GridDataPackageType::
                                            template probeDataPackage<DataType>(pkg_data_addrs, position)
                                      : *pkg_data_addrs[0][0][0];
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
#endif // MESH_WITH_DATA_PACKAGES_3D_HPP