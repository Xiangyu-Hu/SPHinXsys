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
class FinishDataPackages
{
  public:
    explicit FinishDataPackages(MeshWithGridDataPackagesType &mesh_data, Shape &shape)
        : mesh_data_(mesh_data), shape_(shape) {};
    virtual ~FinishDataPackages() {};

    void exec()
    {
        initialize_basic_data_for_a_package.exec();

        near_interface_cell_tagging.exec();
        while (sv_count_modified_.getValue() > 0)
        {
            sv_count_modified_.setValue(0);
            cell_contain_diffusion.exec();
        }

        initialize_cell_neighborhood.exec();
    };

  private:
    MeshWithGridDataPackagesType &mesh_data_;
    Shape &shape_;
    SingularVariable<UnsignedInt> sv_count_modified_{"CountModifiedCell", 1};

    MeshInnerDynamics<execution::ParallelPolicy, InitializeCellNeighborhood> initialize_cell_neighborhood{mesh_data_};
    MeshInnerDynamics<execution::ParallelPolicy, InitializeBasicPackageData> initialize_basic_data_for_a_package{mesh_data_, shape_};
    MeshInnerDynamics<execution::ParallelPolicy, NearInterfaceCellTagging> near_interface_cell_tagging{mesh_data_};
    MeshAllDynamics<execution::ParallelPolicy, CellContainDiffusion> cell_contain_diffusion{mesh_data_, sv_count_modified_};
};

class RepeatTimes
{
  public:
    explicit RepeatTimes() {};
    virtual ~RepeatTimes() {};
    void operator()(UnsignedInt repeat_times) { repeat_times_ = repeat_times; }

  protected:
    UnsignedInt repeat_times_ = 1;
};

template <class ExecutionPolicy>
class CleanInterface : public RepeatTimes, public BaseMeshDynamics, public BaseDynamics<void>
{
  public:
    explicit CleanInterface(MeshWithGridDataPackagesType &mesh_data,
                            NeighborMethod<SingleValued> &neighbor_method,
                            Real refinement_ratio)
        : RepeatTimes(), BaseMeshDynamics(mesh_data), BaseDynamics<void>(),
          neighbor_method_(neighbor_method), refinement_ratio_(refinement_ratio) {};
    virtual ~CleanInterface() {};

    void exec(Real dt = 0.0) override
    {
        for (UnsignedInt k = 0; k != 2 * repeat_times_; ++k)
        {
            for (UnsignedInt l = 0; l != 2; ++l)
            {
                mark_cut_interfaces.exec();
                redistance_interface.exec();
            }

            for (UnsignedInt m = 0; m != 10; ++m)
            {
                reinitialize_level_set.exec();
            }
        }
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    }

  private:
    NeighborMethod<SingleValued> &neighbor_method_;
    Real refinement_ratio_;
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, neighbor_method_};
    MeshInnerDynamics<ExecutionPolicy, MarkCutInterfaces> mark_cut_interfaces{mesh_data_, 0.5 * refinement_ratio_};
    MeshCoreDynamics<ExecutionPolicy, RedistanceInterface> redistance_interface{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, ReinitializeLevelSet> reinitialize_level_set{mesh_data_};
};

template <class ExecutionPolicy>
class CorrectTopology : public BaseMeshDynamics, public BaseDynamics<void>
{
  public:
    explicit CorrectTopology(MeshWithGridDataPackagesType &mesh_data, NeighborMethod<SingleValued> &neighbor_method)
        : BaseMeshDynamics(mesh_data),
          BaseDynamics<void>(),
          neighbor_method_(neighbor_method) {};
    virtual ~CorrectTopology() {};

    void exec(Real dt = 0.0) override
    {
        mark_near_interface.exec();
        while (sv_count_modified_.getValue() > 0)
        {
            sv_count_modified_.setValue(0);
            diffuse_level_set_sign.exec();
        }
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();
    }

  private:
    NeighborMethod<SingleValued> &neighbor_method_;
    SingularVariable<UnsignedInt> sv_count_modified_{"CountModifiedData", 1};
    MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals> update_kernel_integrals{mesh_data_, neighbor_method_};
    MeshInnerDynamics<ExecutionPolicy, MarkNearInterface> mark_near_interface{mesh_data_};
    MeshInnerDynamics<ExecutionPolicy, DiffuseLevelSetSign> diffuse_level_set_sign{mesh_data_, sv_count_modified_};
};
} // namespace SPH
#endif // ALL_MESH_DYNAMICS_H