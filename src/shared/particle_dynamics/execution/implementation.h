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
 * @file 	implementation.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef IMPLEMENTATION_H
#define IMPLEMENTATION_H

#include "base_data_type.h"
#include "base_implementation.h"
#include "execution_policy.h"
#include "ownership.h"

namespace SPH
{
namespace execution
{
template <class ComputingKernelType, class ExecutionPolicy>
inline ComputingKernelType *allocateComputingKernel(const ExecutionPolicy &ex_policy)
{
    return (ComputingKernelType *)malloc(sizeof(ComputingKernelType));
}

template <class ExecutionPolicy, class ComputingKernelType>
inline void copyComputingKernel(const ExecutionPolicy &ex_policy,
                                ComputingKernelType *temp_kernel,
                                ComputingKernelType *computing_kernel)
{
    *computing_kernel = *temp_kernel;
}

template <class ExecutionPolicy, class ComputingKernelType>
inline void freeComputingKernel(const ExecutionPolicy &ex_policy,
                                ComputingKernelType *computing_kernel)
{
    free(computing_kernel);
}

template <class ComputingKernelType>
inline ComputingKernelType *allocateComputingKernelOnDevice();

template <class ComputingKernelType>
inline void copyComputingKernelToDevice(ComputingKernelType *host_kernel,
                                        ComputingKernelType *device_kernel);

template <class ComputingKernelType>
inline void freeComputingKernelOnDevice(ComputingKernelType *device_kernel);

template <class ComputingKernelType, class PolicyType>
inline ComputingKernelType *allocateComputingKernel(const DeviceExecution<PolicyType> &ex_policy)
{
    return allocateComputingKernelOnDevice<ComputingKernelType>();
}

template <class PolicyType, class ComputingKernelType>
inline void copyComputingKernel(const DeviceExecution<PolicyType> &ex_policy,
                                ComputingKernelType *temp_kernel,
                                ComputingKernelType *computing_kernel)
{
    copyComputingKernelToDevice(temp_kernel, computing_kernel);
}

template <class PolicyType, class ComputingKernelType>
inline void freeComputingKernel(const DeviceExecution<PolicyType> &ex_policy,
                                ComputingKernelType *computing_kernel)
{
    freeComputingKernelOnDevice(computing_kernel);
}

template <class ExecutionPolicy, class LocalDynamicsType, class ComputingKernelType>
class Implementation<ExecutionPolicy, LocalDynamicsType, ComputingKernelType>
    : public Implementation<Base>
{
    UniquePtrKeeper<ComputingKernelType> kernel_keeper_;

  public:
    explicit Implementation(LocalDynamicsType &local_dynamics)
        : Implementation<Base>(), local_dynamics_(local_dynamics),
          computing_kernel_(nullptr) {}
    ~Implementation()
    {
        freeComputingKernel(ExecutionPolicy{}, computing_kernel_);
    }

    template <typename... Args>
    ComputingKernelType *getComputingKernel(Args &&...args)
    {
        if (computing_kernel_ == nullptr)
        {
            computing_kernel_ = allocateComputingKernel<ComputingKernelType>(ExecutionPolicy{});
            ComputingKernelType *temp_kernel =
                kernel_keeper_.template createPtr<ComputingKernelType>(
                    ExecutionPolicy{}, this->local_dynamics_, std::forward<Args>(args)...);
            copyComputingKernel(ExecutionPolicy{}, temp_kernel, computing_kernel_);
            this->setUpdated();
        }

        if (!this->isUpdated())
        {
            overwriteComputingKernel(std::forward<Args>(args)...);
        }

        return computing_kernel_;
    }

    template <typename... Args>
    void overwriteComputingKernel(Args &&...args)
    {
        ComputingKernelType *temp_kernel =
            kernel_keeper_.template createPtr<ComputingKernelType>(
                ExecutionPolicy{}, this->local_dynamics_, std::forward<Args>(args)...);
        copyComputingKernel(ExecutionPolicy{}, temp_kernel, computing_kernel_);
        this->setUpdated();
    }

  protected:
    LocalDynamicsType &local_dynamics_;
    ComputingKernelType *computing_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // IMPLEMENTATION_H
