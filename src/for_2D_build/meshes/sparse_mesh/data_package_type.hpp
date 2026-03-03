#ifndef DATA_PACKAGE_TYPE_2D_HPP
#define DATA_PACKAGE_TYPE_2D_HPP

#include "data_package_type.h"

namespace SPH
{
using PackageDataPair = std::pair<UnsignedInt, Array2i>; /**< stores shifted neighbour info: (UnsignedInt)package index, (arrayi)local grid index. */

template <class DataType, int PKG_SIZE>
using PackageBase = PackageBase2d<DataType, PKG_SIZE>;
} // namespace SPH
#endif // DATA_PACKAGE_TYPE_2D_HPP