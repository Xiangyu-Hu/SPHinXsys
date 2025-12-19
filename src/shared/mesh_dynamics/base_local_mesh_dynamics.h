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
#include "spares_mesh_field.hpp"

namespace SPH
{
template <typename DataType>
using PackageVariableData = SparseMeshField<4>::PackageVariableData<DataType>;

template <typename DataType>
using PackageVariable = SparseMeshField<4>::PackageVariable<DataType>;

template <typename DataType>
using CellVariable = DiscreteVariable<DataType>;

template <typename DataType>
using MetaVariable = SparseMeshField<4>::MetaVariable<DataType>;

using PackageVariableAssemble = SparseMeshField<4>::PackageVariableAssemble;
using MetaVariableAssemble = SparseMeshField<4>::MetaVariableAssemble;
using IndexHandler = PackageMesh<4>;

/**
 * @class BaseMeshLocalDynamics
 * @brief The base class for all mesh local particle dynamics.
 */
class BaseMeshLocalDynamics
{
  public:
    explicit BaseMeshLocalDynamics(SparseMeshField<4> &data_mesh, UnsignedInt resolution_level)
        : data_mesh_(data_mesh), index_handler_(data_mesh.getMesh(resolution_level)),
          resolution_level_(resolution_level) {};
    virtual ~BaseMeshLocalDynamics() {};

    static constexpr int pkg_size = SparseMeshField<4>::PackageDataSize();
    static constexpr int pkg_size_minus1 = pkg_size - 1;
    virtual void setupDynamics(Real dt = 0.0) {};  // setup global parameters
    virtual void finishDynamics(Real dt = 0.0) {}; // update global parameters

  protected:
    SparseMeshField<4> &data_mesh_;
    IndexHandler &index_handler_;
    UnsignedInt resolution_level_;
};
} // namespace SPH
#endif // BASE_LOCAL_MESH_DYNAMICS_H
