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
 * @file 	mesh_dynamics.h
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MESH_DYNAMICS_H
#define MESH_DYNAMICS_H

#include "base_dynamics.h"
#include "implementation.h"
#include "mesh_iterators.hpp"
#include "mesh_local_dynamics.hpp"
#include "mesh_with_data_packages.h"
#include "base_configuration_dynamics.h"

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
          all_cells_(mesh_data.AllCells()),
          num_singular_pkgs_(mesh_data.NumSingularPackages()) {};
    virtual ~BaseMeshDynamics() {};

  protected:
    MeshWithGridDataPackagesType &mesh_data_;
    Arrayi all_cells_;
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
        mesh_for(ExecutionPolicy(), MeshRange(Arrayi::Zero(), all_cells_),
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

  public:
    template <typename... Args>
    MeshInnerDynamics(MeshWithGridDataPackagesType &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data), kernel_implementation_(*this){};
    virtual ~MeshInnerDynamics() {};

    template <typename... Args>
    void exec(Args &&...args)
    {
        UnsignedInt num_grid_pkgs = mesh_data_.NumGridPackages();
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        package_for(ExecutionPolicy(), num_singular_pkgs_, num_grid_pkgs,
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
    int *pkg_type_;
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

  public:
    template <typename... Args>
    MeshCoreDynamics(MeshWithGridDataPackagesType &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data),
          pkg_type_(mesh_data.getPackageType().DelegatedData(ExecutionPolicy())),
          kernel_implementation_(*this){};
    virtual ~MeshCoreDynamics() {};

    void exec()
    {
        UnsignedInt num_grid_pkgs = mesh_data_.NumGridPackages();
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        package_for(ExecutionPolicy(), num_singular_pkgs_, num_grid_pkgs,
                    [=](UnsignedInt package_index)
                    {
                        if (pkg_type_[package_index] == 1)
                            update_kernel->update(package_index);
                    });
    };
};

template <class ExecutionPolicy>
class PackageSort : public BaseMeshLocalDynamics, public BaseDynamics<void>
{
    using SortMethodType = typename SortMethod<ExecutionPolicy>::type;

  public:
    explicit PackageSort(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshLocalDynamics(data_mesh), BaseDynamics<void>(),
          occupied_data_pkgs_(data_mesh.getOccupiedDataPackages()),
          dv_sequence_(nullptr), dv_index_permutation_(nullptr),
          update_bkg_mesh_variables_to_sort_(occupied_data_pkgs_.size()),
          update_mesh_variables_to_sort_(occupied_data_pkgs_.size()),
          sort_method_(nullptr) {};
    virtual ~PackageSort() {};
    class ComputingKernel
    {
      public:
        template <class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sequence_(encloser.dv_sequence_.DelegatedData(ex_policy)),
              index_permutation_(encloser.dv_index_permutation_.DelegatedData(ex_policy)){};
        void update(UnsignedInt &pkg_index) {};

      protected:
        UnsignedInt *sequence_, *index_permutation_;
    };

    virtual void exec(Real dt = 0.0) override
    {
        parallel_sort(
            occupied_data_pkgs_.begin(), occupied_data_pkgs_.end(),
            [](const std::pair<UnsignedInt, int> &a, const std::pair<UnsignedInt, int> &b)
            {
                return a.first < b.first;
            });
    };

  private:
    ConcurrentVec<std::pair<UnsignedInt, int>> &occupied_data_pkgs_;
    DiscreteVariable<UnsignedInt> *dv_sequence_;
    DiscreteVariable<UnsignedInt> *dv_index_permutation_;
    OperationOnDataAssemble<BKGMeshVariableAssemble, UpdateSortableVariables<BKGMeshVariable>> update_bkg_mesh_variables_to_sort_;
    OperationOnDataAssemble<MeshVariableAssemble, UpdateSortableVariables<MeshVariable>> update_mesh_variables_to_sort_;
    SortMethodType *sort_method_;
};
} // namespace SPH
#endif // MESH_DYNAMICS_H