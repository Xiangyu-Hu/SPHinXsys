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
namespace SPH
{
namespace execution
{
template <typename... T>
class Implementation;

template <class LocalDynamicsType, class ExecutionPolicy>
class Implementation<ComputingKernelType, ExecutionPolicy>
{
    using ComputingKernel = typename LocalDynamicsType::ComputingKernel;
    UniquePtrKeeper<ComputingKernel> compting_ptr_keeper_;

  public:
    explicit Implementation(LocalDynamicsType &local_dynamics)
        : delegated_kernel_(nullptr) {}

    ComputingKernelType &getDelegatedKernel() const
    {
        if (delegated_kernel_ == nullptr)
        {
            delegated_kernel_ = compting_ptr_keeper_.createPtr<ComputingKernel>(local_dynamics_);
        }
        return *delegated_kernel_;
    }

  private:
    LocalDynamicsType &local_dynamics_;
    ComputingKernelType *delegated_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // EXECUTION_H
