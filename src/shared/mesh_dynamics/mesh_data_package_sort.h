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
 * @file mesh_data_package_sort.h
 * @brief TBD
 * @author Xiangyu Hu
 */

#ifndef MESH_DATA_PACKAGE_SORT_H
#define MESH_DATA_PACKAGE_SORT_H

#include "base_configuration_dynamics.h"
#include "mesh_dynamics_algorithm.h"

namespace SPH
{
template <class ExecutionPolicy>
class PackageSort : public BaseMeshDynamics
{
    using SortMethodType = typename SortMethod<ExecutionPolicy>::type;

  public:
    explicit PackageSort(MeshWithGridDataPackagesType &data_mesh)
        : BaseMeshDynamics(data_mesh),
          ex_policy_(ExecutionPolicy{}),
          sv_num_grid_pkgs_(data_mesh.svNumGridPackages()),
          kernel_implementation_(*this),
          dv_sequence_(data_mesh.registerMetaVariable<UnsignedInt>("Sequence")),
          dv_index_permutation_(data_mesh.registerMetaVariable<UnsignedInt>("IndexPermutation")),
          dv_pkg_1d_cell_index_(&data_mesh.getPackage1DCellIndex()),
          bmv_cell_pkg_index_(&data_mesh.getCellPackageIndex()),
          update_meta_variables_to_sort_(data_mesh.PackageBound()),
          update_mesh_variables_to_sort_(data_mesh.PackageBound()),
          sort_method_(ExecutionPolicy{}, dv_sequence_, dv_index_permutation_) {};
    virtual ~PackageSort() {};

    class UpdateKernel
    {
      public:
        template <class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sequence_(encloser.dv_sequence_->DelegatedData(ex_policy)),
              index_permutation_(encloser.dv_index_permutation_->DelegatedData(ex_policy)),
              pkg_1d_cell_index_(encloser.dv_pkg_1d_cell_index_->DelegatedData(ex_policy)){};
        void update(UnsignedInt &pkg_index)
        {
            sequence_[pkg_index] = pkg_1d_cell_index_[pkg_index];
            index_permutation_[pkg_index] = pkg_index;
        };

      protected:
        UnsignedInt *sequence_, *index_permutation_, *pkg_1d_cell_index_;
    };

    void exec(Real dt = 0.0)
    {
        UnsignedInt num_grid_pkgs = sv_num_grid_pkgs_.getValue();
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        package_for(ex_policy_, 0, num_grid_pkgs,
                    [=](UnsignedInt package_index)
                    {
                        update_kernel->update(package_index);
                    });

        UnsignedInt sortable_size = num_grid_pkgs - num_singular_pkgs_;
        sort_method_.sort(ex_policy_, sortable_size, num_singular_pkgs_);
        update_meta_variables_to_sort_(
            mesh_data_.getEvolvingMetaVariables(),
            ex_policy_, num_grid_pkgs, dv_index_permutation_);
        update_mesh_variables_to_sort_(
            mesh_data_.getEvolvingMeshVariables(),
            ex_policy_, num_grid_pkgs, dv_index_permutation_);

        UnsignedInt *pkg_1d_cell_index = dv_pkg_1d_cell_index_->DelegatedData(ex_policy_);
        UnsignedInt *cell_pkg_index = bmv_cell_pkg_index_->DelegatedData(ex_policy_);
        package_for(ex_policy_, num_singular_pkgs_, num_grid_pkgs,
                    [=](UnsignedInt package_index)
                    {
                        UnsignedInt sort_index = pkg_1d_cell_index[package_index];
                        cell_pkg_index[sort_index] = package_index;
                    });
    };

  private:
    ExecutionPolicy ex_policy_;
    SingularVariable<UnsignedInt> &sv_num_grid_pkgs_;
    using KernelImplementation = Implementation<ExecutionPolicy, PackageSort<ExecutionPolicy>, UpdateKernel>;
    KernelImplementation kernel_implementation_;
    DiscreteVariable<UnsignedInt> *dv_sequence_;
    DiscreteVariable<UnsignedInt> *dv_index_permutation_;
    MetaVariable<UnsignedInt> *dv_pkg_1d_cell_index_;
    BKGMeshVariable<UnsignedInt> *bmv_cell_pkg_index_;
    OperationOnDataAssemble<MetaVariableAssemble, UpdateSortableVariables<MetaVariable>> update_meta_variables_to_sort_;
    OperationOnDataAssemble<MeshVariableAssemble, UpdateSortableVariables<MeshVariable>> update_mesh_variables_to_sort_;
    SortMethodType sort_method_;
};
} // namespace SPH
#endif // MESH_DATA_PACKAGE_SORT_H