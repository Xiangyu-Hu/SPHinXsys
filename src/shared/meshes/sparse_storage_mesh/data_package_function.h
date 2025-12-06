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
 * @file data_package_function.h
 * @brief This is for the base functions for mesh iterator.
 * There are two types of functions: one is for static ranged
 * which are defined by template parameters,
 * the other for dynamics ranges which are input parameters.
 * @author  Xiangyu Hu
 */

#ifndef DATA_PACKAGE_FUNCTIONS_H
#define DATA_PACKAGE_FUNCTIONS_H

#include "mesh_with_data_packages.h"

namespace SPH
{
template <int PKG_SIZE>
DataPackagePair NeighbourIndexShift(const Arrayi &shift_index, const CellNeighborhood &neighbour);

template <int PKG_SIZE>
DataPackagePair GeneralNeighbourIndexShift(
    UnsignedInt package_index, CellNeighborhood *neighbour, const Arrayi &shift_index);

template <typename DataType, int PKG_SIZE>
DataType CornerAverage(PackageData<DataType, PKG_SIZE> *pkg_data, Arrayi addrs_index,
                       Arrayi corner_direction, const CellNeighborhood &neighborhood, DataType zero);

template <typename DataType, int PKG_SIZE, typename FunctionByIndex>
void assignByDataIndex(PackageData<DataType, PKG_SIZE> &pkg_data, const FunctionByIndex &function_by_index);

template <int PKG_SIZE, typename RegularizeFunction>
Vec2d regularizedCentralDifference(PackageData<Real, PKG_SIZE> *input, const CellNeighborhood2d &neighborhood,
                                   const Array2i &data_index, const RegularizeFunction &regularize_function);

template <int PKG_SIZE, typename RegularizeFunction>
Vec3d regularizedCentralDifference(PackageData<Real, PKG_SIZE> *input, const CellNeighborhood3d &neighborhood,
                                   const Array3i &data_index, const RegularizeFunction &regularize_function);
template <typename DataType, int PKG_SIZE>
class ProbeMesh
{
    using IndexHandler = typename MeshWithGridDataPackages<PKG_SIZE>::IndexHandler;

  public:
    template <class ExecutionPolicy>
    ProbeMesh(const ExecutionPolicy &ex_policy, MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
              const std::string variable_name);
    DataType operator()(const Vecd &position);

  protected:
    PackageData<DataType, PKG_SIZE> *pkg_data_;
    IndexHandler index_handler_;
    UnsignedInt *cell_pkg_index_;
    CellNeighborhood *cell_neighborhood_;
    /** probe by applying bi and tri-linear interpolation within the package. */
    DataType probeDataPackage(UnsignedInt package_index, const Array2i &cell_index, const Vec2d &position);
    DataType probeDataPackage(UnsignedInt package_index, const Array3i &cell_index, const Vec3d &position);
};
} // namespace SPH
#endif // DATA_PACKAGE_FUNCTIONS_H