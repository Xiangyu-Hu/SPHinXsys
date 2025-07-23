#ifndef GRID_DATA_PACKAGE_TYPE_2D_HPP
#define GRID_DATA_PACKAGE_TYPE_2D_HPP

#include "grid_data_package_type.h"

namespace SPH
{
using CellNeighborhood = CellNeighborhood2d;
using NeighbourIndex = std::pair<size_t, Array2i>; /**< stores shifted neighbour info: (size_t)package index, (arrayi)local grid index. */

template <class DataType, size_t PKG_SIZE>
using PackageDataMatrix = PackageDataMatrix2d<DataType, PKG_SIZE>;
} // namespace SPH
#endif // GRID_DATA_PACKAGE_TYPE_2D_HPP