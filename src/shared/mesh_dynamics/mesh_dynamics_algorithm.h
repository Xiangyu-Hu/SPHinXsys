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
#include "mesh_with_data_packages.h"

namespace SPH
{
/**
 * @class BaseMeshDynamics
 * @brief The base class for all mesh dynamics
 * This class contains only the interface functions available
 * for all dynamics. An specific implementation should be realized.
 */
class BaseMeshDynamics
{
  public:
    BaseMeshDynamics(MeshWithGridDataPackagesType &mesh_data)
        : mesh_data_(mesh_data),
          index_handler_(mesh_data_.getIndexHandler()),
          num_singular_pkgs_(mesh_data.NumSingularPackages()) {};
    virtual ~BaseMeshDynamics() {};

  protected:
    MeshWithGridDataPackagesType &mesh_data_;
    IndexHandler &index_handler_;
    UnsignedInt num_singular_pkgs_;
};

/**
 * @class MeshAllDynamics
 * @brief Mesh dynamics for all cell on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshAllDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    MeshAllDynamics(MeshWithGridDataPackagesType &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data), kernel_implementation_(*this){};
    virtual ~MeshAllDynamics() {};

    void exec()
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        mesh_for(ExecutionPolicy(), MeshRange(Arrayi::Zero(), index_handler_.AllCells()),
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
class MeshInnerDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
    SingularVariable<UnsignedInt> &sv_num_grid_pkgs_;

  public:
    template <typename... Args>
    MeshInnerDynamics(MeshWithGridDataPackagesType &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data), kernel_implementation_(*this),
          sv_num_grid_pkgs_(mesh_data.svNumGridPackages()){};
    virtual ~MeshInnerDynamics() {};

    template <typename... Args>
    void exec(Args &&...args)
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        package_for(ExecutionPolicy(), num_singular_pkgs_, sv_num_grid_pkgs_.getValue(),
                    [=](UnsignedInt package_index)
                    {
                        update_kernel->update(package_index, args...);
                    });
    };
};

/**
 * @class MeshCoreDynamics
 * @brief Mesh dynamics for only core cells on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshCoreDynamics : public LocalDynamicsType, public BaseMeshDynamics
{
    MetaVariable<int> &dv_pkg_type_;
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
    SingularVariable<UnsignedInt> &sv_num_grid_pkgs_;

  public:
    template <typename... Args>
    MeshCoreDynamics(MeshWithGridDataPackagesType &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data), dv_pkg_type_(mesh_data.getPackageType()),
          kernel_implementation_(*this),
          sv_num_grid_pkgs_(mesh_data.svNumGridPackages()){};
    virtual ~MeshCoreDynamics() {};

    void exec()
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        int *pkg_type = dv_pkg_type_.DelegatedData(ExecutionPolicy());
        package_for(ExecutionPolicy(), num_singular_pkgs_, sv_num_grid_pkgs_.getValue(),
                    [=](UnsignedInt package_index)
                    {
                        if (pkg_type[package_index] == 1)
                            update_kernel->update(package_index);
                    });
    };
};
} // namespace SPH
#endif // MESH_DYNAMICS_ALGORITHM_H