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
 * @file 	all_mesh_dynamics.h
 * @brief 	This is for the base classes of particle dynamics, which describe the
 * 			interaction between particles. These interactions are used to define
 *			differential operators for surface forces or fluxes in continuum mechanics
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALL_MESH_DYNAMICS_H
#define ALL_MESH_DYNAMICS_H

#include "mesh_dynamics.h"
#include "mesh_local_dynamics.hpp"

namespace SPH
{
class BaseExecDynamics
{
  public:
    BaseExecDynamics() {};
    virtual ~BaseExecDynamics() {};

    virtual void exec(Real small_shift_factor) = 0;
};
class FinishDataPackages
{
  public:
    explicit FinishDataPackages(MeshWithGridDataPackagesType &mesh_data, Shape &shape)
        : mesh_data_(mesh_data), shape_(shape),
          far_field_distance(mesh_data.GridSpacing() * (Real)mesh_data.BufferWidth()) {};
    virtual ~FinishDataPackages() {};

    void exec()
    {
        tag_a_cell_is_inner_package.exec();

        mesh_data_.organizeOccupiedPackages();
        initialize_index_mesh.exec();
        initialize_cell_neighborhood.exec();
        mesh_data_.resizeMeshVariableData();

        initialize_data_for_singular_package.update(0, -far_field_distance);
        initialize_data_for_singular_package.update(1, far_field_distance);

        initialize_basic_data_for_a_package.exec();
    };

  private:
    MeshWithGridDataPackagesType &mesh_data_;
    Shape &shape_;
    Real far_field_distance;

    InitializeDataForSingularPackage initialize_data_for_singular_package{mesh_data_};
    MeshAllDynamics<execution::ParallelPolicy, InnerCellITagging> tag_a_cell_is_inner_package{mesh_data_};
    MeshInnerDynamics<execution::ParallelPolicy, InitializeIndexMesh> initialize_index_mesh{mesh_data_};
    MeshInnerDynamics<execution::ParallelPolicy, InitializeCellNeighborhood> initialize_cell_neighborhood{mesh_data_};
    MeshInnerDynamics<execution::ParallelPolicy, InitializeBasicDataForAPackage> initialize_basic_data_for_a_package{mesh_data_, shape_};
};

template <class ExecutionPolicy>
class CleanInterface : public BaseMeshDynamics, public BaseExecDynamics
{
  public:
    explicit CleanInterface(MeshWithGridDataPackagesType &mesh_data, KernelTabulatedCK *kernel, Real global_h_ratio)
        : BaseMeshDynamics(mesh_data),
          BaseExecDynamics(),
          kernel_(kernel),
          global_h_ratio_(global_h_ratio) {};
    virtual ~CleanInterface() {};

    void exec(Real small_shift_factor) override
    {
        mark_near_interface.exec(small_shift_factor);
        redistance_interface.exec();
        reinitialize_level_set.exec();
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    }

  private:
    KernelTabulatedCK *kernel_;
    Real global_h_ratio_;
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, kernel_, global_h_ratio_};
    MeshInnerDynamics<ExecutionPolicy, MarkNearInterface> mark_near_interface{mesh_data_};
    MeshCoreDynamics<ExecutionPolicy, RedistanceInterface> redistance_interface{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, ReinitializeLevelSet> reinitialize_level_set{mesh_data_};
};

template <class ExecutionPolicy>
class CorrectTopology : public BaseMeshDynamics, public BaseExecDynamics
{
  public:
    explicit CorrectTopology(MeshWithGridDataPackagesType &mesh_data, KernelTabulatedCK *kernel, Real global_h_ratio)
        : BaseMeshDynamics(mesh_data),
          BaseExecDynamics(),
          kernel_(kernel),
          global_h_ratio_(global_h_ratio) {};
    virtual ~CorrectTopology() {};

    void exec(Real small_shift_factor) override
    {
        mark_near_interface.exec(small_shift_factor);
        for (size_t i = 0; i != 10; ++i)
            diffuse_level_set_sign.exec();
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    }

  private:
    KernelTabulatedCK *kernel_;
    Real global_h_ratio_;
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, kernel_, global_h_ratio_};
    MeshInnerDynamics<ExecutionPolicy, MarkNearInterface> mark_near_interface{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, DiffuseLevelSetSign> diffuse_level_set_sign{mesh_data_};
};
} // namespace SPH
#endif // ALL_MESH_DYNAMICS_H