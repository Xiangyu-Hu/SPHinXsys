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
 * @file    base_local_mesh_dynamics.h
 * @brief   TBD.
 * @author  Xiangyu Hu
 */

#ifndef BASE_LOCAL_MESH_DYNAMICS_H
#define BASE_LOCAL_MESH_DYNAMICS_H

#include "base_dynamics.h"
#include "base_implementation.h"
#include "data_package_function.hpp"
#include "mesh_with_data_packages.hpp"

namespace SPH
{
using MeshWithGridDataPackagesType = MeshWithGridDataPackages<4>;

template <typename DataType>
using MeshVariableData = MeshWithGridDataPackagesType::MeshVariableData<DataType>;

template <typename DataType>
using MeshVariable = MeshWithGridDataPackagesType::MeshVariable<DataType>;

template <typename DataType>
using BKGMeshVariable = MeshWithGridDataPackagesType::BKGMeshVariable<DataType>;

template <typename DataType>
using MetaVariable = MeshWithGridDataPackagesType::MetaVariable<DataType>;

using MeshVariableAssemble = MeshWithGridDataPackagesType::MeshVariableAssemble;
using BKGMeshVariableAssemble = MeshWithGridDataPackagesType::BKGMeshVariableAssemble;
using MetaVariableAssemble = MeshWithGridDataPackagesType::MetaVariableAssemble;
using IndexHandler = MeshWithGridDataPackagesType::IndexHandler;

/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(MeshWithGridDataPackagesType &data_mesh)
        : data_mesh_(data_mesh), index_handler_(data_mesh.getIndexHandler()) {};
    virtual ~BaseMeshLocalDynamics() {};

    MeshWithGridDataPackagesType &data_mesh_;
    IndexHandler &index_handler_;
    static constexpr int pkg_size = MeshWithGridDataPackagesType::DataPackageSize();
    static constexpr int pkg_size_minus1 = pkg_size - 1;
};
} // namespace SPH
#endif // BASE_LOCAL_MESH_DYNAMICS_H
