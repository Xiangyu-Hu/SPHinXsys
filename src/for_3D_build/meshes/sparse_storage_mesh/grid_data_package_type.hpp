#ifndef GRID_DATA_PACKAGE_TYPE_3D_HPP
#define GRID_DATA_PACKAGE_TYPE_3D_HPP

#include "grid_data_package_type.h"

namespace SPH
{
using CellNeighborhood = CellNeighborhood3d;
using PackageGridPair = std::pair<UnsignedInt, Array3i>; /**< stores shifted neighbour info: (UnsignedInt)package index, (arrayi)local grid index. */

template <class DataType, UnsignedInt PKG_SIZE>
using PackageDataMatrix = PackageDataMatrix3d<DataType, PKG_SIZE>;
} // namespace SPH
#endif // GRID_DATA_PACKAGE_TYPE_3D_HPP