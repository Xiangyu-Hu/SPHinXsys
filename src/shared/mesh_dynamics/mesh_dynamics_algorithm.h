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
 * @file 	mesh_dynamics_algorithm.h
 * @brief TBD
 * @author Xiangyu Hu
 */

#ifndef MESH_DYNAMICS_ALGORITHM_H
#define MESH_DYNAMICS_ALGORITHM_H

#include "base_configuration_dynamics.h"
#include "base_dynamics.h"
#include "base_local_mesh_dynamics.h"
#include "implementation.h"
#include "mesh_iterators.hpp"
#include "spares_mesh_field.h"

namespace SPH
{
/**
 * @class MeshAllDynamics
 * @brief Mesh dynamics for all cell on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshAllDynamics : public LocalDynamicsType, public BaseDynamics<void>
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    MeshAllDynamics(
        SparseMeshField<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseDynamics<void>(), kernel_implementation_(*this){};
    virtual ~MeshAllDynamics() {};

    void exec(Real dt = 0.0) override
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        mesh_for(ExecutionPolicy(), MeshRange(Arrayi::Zero(), this->index_handler_.AllCells()),
                 [&](Arrayi cell_index)
                 {
                     update_kernel->update(cell_index);
                 });
    };
};

/**
 * @class MeshInnerDynamics
 * @brief Mesh dynamics for only inner cells on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshInnerDynamics : public LocalDynamicsType, public BaseDynamics<void>
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
    StdVec<UnsignedInt> &num_pkgs_offsets_;

  public:
    template <typename... Args>
    MeshInnerDynamics(
        SparseMeshField<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseDynamics<void>(), kernel_implementation_(*this),
          num_pkgs_offsets_(mesh_data.getNumPackageOffsets()){};
    virtual ~MeshInnerDynamics() {};

    void exec(Real dt = 0.0) override
    {
        LocalDynamicsType::setupDynamics(dt);
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        UnsignedInt start_pkg_index = num_pkgs_offsets_[this->resolution_level_];
        UnsignedInt end_pkg_index = num_pkgs_offsets_[this->resolution_level_ + 1];
        package_for(ExecutionPolicy(), start_pkg_index, end_pkg_index,
                    [=](UnsignedInt package_index)
                    {
                        update_kernel->update(package_index);
                    });
        LocalDynamicsType::finishDynamics(dt);
    };
};

/**
 * @class MeshCoreDynamics
 * @brief Mesh dynamics for only core cells on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshCoreDynamics : public LocalDynamicsType, public BaseDynamics<void>
{
    MetaVariable<int> &dv_pkg_type_;
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
    StdVec<UnsignedInt> &num_pkgs_offsets_;

  public:
    template <typename... Args>
    MeshCoreDynamics(
        SparseMeshField<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseDynamics<void>(), dv_pkg_type_(mesh_data.getPackageType()),
          kernel_implementation_(*this),
          num_pkgs_offsets_(mesh_data.getNumPackageOffsets()){};
    virtual ~MeshCoreDynamics() {};

    void exec(Real dt = 0.0) override
    {
        LocalDynamicsType::setupDynamics(dt);
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        int *pkg_type = dv_pkg_type_.DelegatedData(ExecutionPolicy());
        UnsignedInt start_pkg_index = num_pkgs_offsets_[this->resolution_level_];
        UnsignedInt end_pkg_index = num_pkgs_offsets_[this->resolution_level_ + 1];
        package_for(ExecutionPolicy(), start_pkg_index, end_pkg_index,
                    [=](UnsignedInt package_index)
                    {
                        if (pkg_type[package_index] == 1)
                            update_kernel->update(package_index);
                    });
        LocalDynamicsType::finishDynamics(dt);
    };
};
} // namespace SPH
#endif // MESH_DYNAMICS_ALGORITHM_H