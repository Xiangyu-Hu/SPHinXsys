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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
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

#include "mesh_local_dynamics.hpp"
#include "mesh_with_data_packages.h"
#include "mesh_iterators.hpp"
#include "execution.h"

#include <functional>

using namespace std::placeholders;

namespace SPH
{

/***********************************************
 *       Iterators for Only Occupied Cells      *
 ***********************************************/
template <typename FunctionOnData>
void package_parallel_for(const execution::SequencedPolicy &seq,
                          size_t num_grid_pkgs, const FunctionOnData &function)
{
    for (size_t i = 2; i != num_grid_pkgs; ++i)
        function(i);
}
template <typename FunctionOnData>
void package_parallel_for(const execution::ParallelPolicy &par,
                          size_t num_grid_pkgs, const FunctionOnData &function)
{
    parallel_for(IndexRange(2, num_grid_pkgs),
                [&](const IndexRange &r)
                {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        function(i);
                    }
                },
                ap);
}
template <typename FunctionOnData>
void package_parallel_for(const execution::ParallelDevicePolicy &par_device,
                          size_t num_grid_pkgs, const FunctionOnData &function);
/**
 * @class BaseMeshDynamics
 * @brief The base class for all mesh dynamics
 * This class contains only the interface functions available
 * for all dynamics. An specific implementation should be realized.
 */
class BaseMeshDynamics
{
  public:
    BaseMeshDynamics(MeshWithGridDataPackages<4> &mesh_data)
        : mesh_data_(mesh_data),
          all_cells_(mesh_data.AllCells()),
          num_grid_pkgs_(mesh_data.num_grid_pkgs_){};
    virtual ~BaseMeshDynamics(){};

  protected:
    MeshWithGridDataPackages<4> &mesh_data_;
    Arrayi all_cells_;
    size_t &num_grid_pkgs_;

    /***********************************************
     *         Iterators for All Mesh Cells        *
     ***********************************************/
    template <typename FunctionOnData>
    void grid_parallel_for(const execution::SequencedPolicy &seq, const FunctionOnData &function)
    {
        mesh_for(MeshRange(Arrayi::Zero(), all_cells_),
                          [&](Arrayi cell_index)
                          {
                              function(cell_index);
                          });
    }
    template <typename FunctionOnData>
    void grid_parallel_for(const execution::ParallelPolicy &par, const FunctionOnData &function)
    {
        mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                          [&](Arrayi cell_index)
                          {
                              function(cell_index);
                          });
    }

    /***********************************************
     *       Iterators for Only Occupied Cells      *
     ***********************************************/
    // template <typename FunctionOnData>
    // void package_parallel_for(const execution::SequencedPolicy &seq, const FunctionOnData &function)
    // {
    //     for (size_t i = 2; i != num_grid_pkgs_; ++i)
    //         function(i);
    // }
    // template <typename FunctionOnData>
    // void package_parallel_for(const execution::ParallelPolicy &par, const FunctionOnData &function)
    // {
    //     parallel_for(IndexRange(2, num_grid_pkgs_),
    //                 [&](const IndexRange &r)
    //                 {
    //                     for (size_t i = r.begin(); i != r.end(); ++i)
    //                     {
    //                         function(i);
    //                     }
    //                 },
    //                 ap);
    // }

    // #if SPHINXSYS_USE_SYCL
    // template <typename FunctionOnData>
    // void BaseMeshDynamics::package_parallel_for(const execution::ParallelDevicePolicy &par_device, const FunctionOnData &function)
    // {
    //     auto &sycl_queue = execution_instance.getQueue();
    //     sycl_queue.submit([&](sycl::handler &cgh)
    //     {
    //         const size_t num_grid_pkgs = this->num_grid_pkgs_;
    //         cgh.parallel_for(execution_instance.getUniformNdRange(num_grid_pkgs), [=](sycl::nd_item<1> index)
    //         {
    //             if(index.get_global_id(0) + 2 < num_grid_pkgs)
    //                 function(index.get_global_id(0) + 2); 
    //         });
    //     }).wait_and_throw();
    // }
    // #endif

    // template <typename FunctionOnData>
    // void package_parallel_for(const execution::ParallelDevicePolicy &par_device, const FunctionOnData &function);
    // {
    //     auto &sycl_queue = execution_instance.getQueue();
    //     sycl_queue.submit([&](sycl::handler &cgh)
    //                       { cgh.parallel_for(execution_instance.getUniformNdRange(num_grid_pkgs_ - 2), [=](sycl::nd_item<1> index)
    //                                         {
    //                                 if(index.get_global_id(0) + 2< num_grid_pkgs_)
    //                                     local_dynamics_function(index.get_global_id(0) + 2); }); })
    //         .wait_and_throw();
    // }
};

/**
 * @class MeshAllDynamics
 * @brief Mesh dynamics for all cell on the mesh (SYCL version)
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshAllDynamicsCK : public LocalDynamicsType, public BaseMeshDynamics
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
  public:
    template <typename... Args>
    MeshAllDynamicsCK(MeshWithGridDataPackages<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data), kernel_implementation_(*this){};
    virtual ~MeshAllDynamicsCK(){};

    void exec()
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        grid_parallel_for(ExecutionPolicy(),
            [&](Arrayi cell_index)
            {
              update_kernel->update(cell_index);
            }
        );
    };
};

/**
 * @class MeshInnerDynamics
 * @brief Mesh dynamics for only inner cells on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshInnerDynamicsCK : public LocalDynamicsType, public BaseMeshDynamics
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;
  public:
    template <typename... Args>
    MeshInnerDynamicsCK(MeshWithGridDataPackages<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data), kernel_implementation_(*this){};
    virtual ~MeshInnerDynamicsCK(){};

    template <typename... Args>
    void exec(Args &&...args)
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        package_parallel_for(ExecutionPolicy(), num_grid_pkgs_,
            [=](size_t package_index)
            {
              update_kernel->update(package_index, args...);
            }
        );
    };
};

/**
 * @class MeshCoreDynamics
 * @brief Mesh dynamics for only core cells on the mesh
 */
template <class ExecutionPolicy, class LocalDynamicsType>
class MeshCoreDynamicsCK : public LocalDynamicsType, public BaseMeshDynamics
{
    using UpdateKernel = typename LocalDynamicsType::UpdateKernel;
    using KernelImplementation = Implementation<ExecutionPolicy, LocalDynamicsType, UpdateKernel>;
    KernelImplementation kernel_implementation_;

    std::pair<Arrayi, int> *meta_data_cell_;
  public:
    template <typename... Args>
    MeshCoreDynamicsCK(MeshWithGridDataPackages<4> &mesh_data, Args &&...args)
        : LocalDynamicsType(mesh_data, std::forward<Args>(args)...),
          BaseMeshDynamics(mesh_data),
          meta_data_cell_(mesh_data.meta_data_cell_.DelegatedDataField(ExecutionPolicy())),
          kernel_implementation_(*this){};
    virtual ~MeshCoreDynamicsCK(){};

    void exec()
    {
        UpdateKernel *update_kernel = kernel_implementation_.getComputingKernel();
        std::pair<SPH::Arrayi, int> *meta_data_cell = meta_data_cell_;
        package_parallel_for(ExecutionPolicy(), num_grid_pkgs_,
            [=](size_t package_index)
            {
              if (meta_data_cell[package_index].second == 1)
                  update_kernel->update(package_index);
            }
        );
    };
};



} // namespace SPH
#endif // MESH_DYNAMICS_H