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
 * @file 	execution.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef EXECUTION_H
#define EXECUTION_H

#include "execution_policy.h"

#include "base_data_type.h"
#include "ownership.h"

namespace SPH
{
namespace execution
{
template <typename... T>
class Implementation;

template <class LocalDynamicsType, class ExecutionPolicy>
class Implementation<LocalDynamicsType, ExecutionPolicy>
{
    using ComputingKernel = typename LocalDynamicsType::
        template ComputingKernel<ExecutionPolicy>;
    UniquePtrKeeper<ComputingKernel> kernel_ptr_keeper_;

  public:
    explicit Implementation(LocalDynamicsType &local_dynamics)
        : local_dynamics_(local_dynamics), computing_kernel_(nullptr) {}
    ~Implementation()
    {
        delete computing_kernel_;
    }

    ComputingKernel *getComputingKernel()
    {
        if (computing_kernel_ == nullptr)
        {
            computing_kernel_ = new ComputingKernel;
            ComputingKernel *temp_kernel =
                kernel_ptr_keeper_.template createPtr<ComputingKernel>(ExecutionPolicy{}, local_dynamics_);
            *computing_kernel_ = *temp_kernel;
        }
        return computing_kernel_;
    }

    ComputingKernel *overwriteComputingKernel()
    {
        ComputingKernel *temp_kernel =
            kernel_ptr_keeper_.template createPtr<ComputingKernel>(ExecutionPolicy{}, local_dynamics_);
        *computing_kernel_ = *temp_kernel;
        return computing_kernel_;
    }

  private:
    LocalDynamicsType &local_dynamics_;
    ComputingKernel *computing_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // EXECUTION_H
