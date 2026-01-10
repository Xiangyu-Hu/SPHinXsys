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
 * @file data_package_utility.h
 * @brief TBD.
 * @author  Xiangyu Hu
 */

#ifndef DATA_PACKAGE_UTILITY_H
#define DATA_PACKAGE_UTILITY_H

#include "data_package_type.hpp"

namespace SPH
{
template <class DataType, int PKG_SIZE>
class PackageData : public PackageBase<DataType, PKG_SIZE>
{
  public:
    static PackageData Constant(const DataType &value)
    {
        PackageData pkg_data{};
        mesh_for_each(Arrayi::Zero(), Arrayi::Constant(PKG_SIZE),
                      [&](const Arrayi &data_index)
                      {
                          pkg_data(data_index) = value;
                      });
        return pkg_data;
    };

    template <typename FunctionByIndex>
    void assignByIndex(const FunctionByIndex &function_by_index);
};

class CellNeighborhood : public PackageData<UnsignedInt, 3>
{
  public:
    CellNeighborhood() = default;
    CellNeighborhood(const PackageData<UnsignedInt, 3> &base_obj)
        : PackageData<UnsignedInt, 3>(base_obj) {};

    template <int PKG_SIZE>
    PackageDataPair IndexShift(const Arrayi &shift_index);

    template <typename DataType, int PKG_SIZE>
    DataType DataFromIndex(PackageData<DataType, PKG_SIZE> *pkg_data, const Arrayi &data_index);

    template <typename DataType, int PKG_SIZE>
    DataType CornerAverage(PackageData<DataType, PKG_SIZE> *pkg_data, Arrayi addrs_index,
                           Arrayi corner_direction, DataType zero);

    template <int PKG_SIZE, typename RegularizeFunction>
    Vec2d regularizedCentralDifference(PackageData<Real, PKG_SIZE> *input, const Array2i &data_index,
                                       const RegularizeFunction &regularize_function);
    template <int PKG_SIZE, typename RegularizeFunction>
    Vec3d regularizedCentralDifference(PackageData<Real, PKG_SIZE> *input, const Array3i &data_index,
                                       const RegularizeFunction &regularize_function);
};

template <int PKG_SIZE>
PackageDataPair GeneralNeighbourIndexShift(
    UnsignedInt package_index, CellNeighborhood *neighbour, const Arrayi &shift_index);
} // namespace SPH
#endif // DATA_PACKAGE_UTILITY_H