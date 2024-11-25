/**
 * @file    mesh_with_data_packages.hpp
 * @brief   Implementation for 2d builds.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_WITH_DATA_PACKAGES_2D_HPP
#define MESH_WITH_DATA_PACKAGES_2D_HPP

#include "mesh_with_data_packages.h"
#include "mesh_iterators.hpp"

namespace SPH
{
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::allocateIndexDataMatrix()
{
    Allocate2dArray(index_data_mesh_, all_cells_);
}
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::deleteIndexDataMatrix()
{
    Delete2dArray(index_data_mesh_, all_cells_);
}
//=================================================================================================//
template <int PKG_SIZE>
void MeshWithGridDataPackages<PKG_SIZE>::
    assignDataPackageIndex(const Arrayi &cell_index, const size_t package_index)
{
    index_data_mesh_[cell_index[0]][cell_index[1]] = package_index;
}
//=================================================================================================//
template <int PKG_SIZE>
size_t MeshWithGridDataPackages<PKG_SIZE>::
    PackageIndexFromCellIndex(const Arrayi &cell_index)
{
    return index_data_mesh_[cell_index[0]][cell_index[1]];
}
//=================================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::
    isSingularDataPackage(const Arrayi &cell_index)
{
    return index_data_mesh_[cell_index[0]][cell_index[1]] < 2;
}
//=================================================================================================//
template <int PKG_SIZE>
bool MeshWithGridDataPackages<PKG_SIZE>::
    isInnerDataPackage(const Arrayi &cell_index)
{
    return index_data_mesh_[cell_index[0]][cell_index[1]] > 1;
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//
#endif // MESH_WITH_DATA_PACKAGES_2D_HPP