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
 * @file grid_data_package_function.h
 * @brief This is for the base functions for mesh iterator.
 * There are two types of functions: one is for static ranged
 * which are defined by template parameters,
 * the other for dynamics ranges which are input parameters.
 * @author  Xiangyu Hu
 */

#ifndef GRID_DATA_PACKAGE_FUNCTIONS_H
#define GRID_DATA_PACKAGE_FUNCTIONS_H

#include "mesh_with_data_packages.h"

namespace SPH
{
template <int PKG_SIZE>
NeighbourIndex NeighbourIndexShift(const Arrayi shift_index, const CellNeighborhood &neighbour);

template <typename DataType, size_t PKG_SIZE>
DataType CornerAverage(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data, Arrayi addrs_index,
                       Arrayi corner_direction, const CellNeighborhood &neighborhood, DataType zero);

template <typename DataType, size_t PKG_SIZE>
DataType DataValueFromGlobalIndex(PackageDataMatrix<DataType, PKG_SIZE> *pkg_data,
                                  const Arrayi &global_grid_index,
                                  MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
                                  size_t *cell_package_index);

template <typename DataType, size_t PKG_SIZE>
class ProbeMesh
{
    using IndexHandler = typename MeshWithGridDataPackages<PKG_SIZE>::IndexHandler;

  public:
    template <class ExecutionPolicy>
    ProbeMesh(const ExecutionPolicy &ex_policy, MeshWithGridDataPackages<PKG_SIZE> *data_mesh,
              const std::string variable_name)
        : pkg_data_(data_mesh->template getMeshVariable<DataType>(variable_name)->DelegatedData(ex_policy)),
          index_handler_(data_mesh->index_handler_.DelegatedData(ex_policy)),
          cell_package_index_(data_mesh->cell_package_index_.DelegatedData(ex_policy)),
          cell_neighborhood_(data_mesh->cell_neighborhood_.DelegatedData(ex_policy)){};

    DataType operator()(const Vecd &position);

  protected:
    PackageDataMatrix<DataType, PKG_SIZE> *pkg_data_;
    IndexHandler *index_handler_;
    size_t *cell_package_index_;
    CellNeighborhood *cell_neighborhood_;
    /** probe by applying bi and tri-linear interpolation within the package. */
    DataType probeDataPackage(size_t package_index, const Arrayi &cell_index, const Vecd &position);
};
} // namespace SPH
#endif // GRID_DATA_PACKAGE_FUNCTIONS_H