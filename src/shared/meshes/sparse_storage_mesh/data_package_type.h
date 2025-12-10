/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file data_package_type.h
 * @brief TBD.
 * @author  Xiangyu Hu
 */

#ifndef DATA_PACKAGE_TYPE_H
#define DATA_PACKAGE_TYPE_H

#include "base_data_type_package.h"

namespace SPH
{
template <int PKG_SIZE>
inline constexpr std::size_t as_size_t_v = static_cast<std::size_t>(PKG_SIZE);

template <class DataType, int PKG_SIZE>
class PackageData2d
    : public std::array<std::array<DataType, as_size_t_v<PKG_SIZE>>, as_size_t_v<PKG_SIZE>>
{
  public:
    DataType operator()(const Array2i &index) const { return (*this)[index[0]][index[1]]; }
    DataType &operator()(const Array2i &index) { return (*this)[index[0]][index[1]]; }
};

template <class DataType, int PKG_SIZE>
class PackageData3d
    : public std::array<std::array<std::array<DataType, as_size_t_v<PKG_SIZE>>, as_size_t_v<PKG_SIZE>>, as_size_t_v<PKG_SIZE>>
{
  public:
    DataType operator()(const Array3i &index) const { return (*this)[index[0]][index[1]][index[2]]; }
    DataType &operator()(const Array3i &index) { return (*this)[index[0]][index[1]][index[2]]; }
};

using CellNeighborhood2d = PackageData2d<UnsignedInt, 3>;
using CellNeighborhood3d = PackageData3d<UnsignedInt, 3>;
} // namespace SPH
#endif // DATA_PACKAGE_TYPE_H