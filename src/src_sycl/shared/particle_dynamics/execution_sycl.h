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
 * @file 	execution_sycl.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef EXECUTION_SYCL_H
#define EXECUTION_SYCL_H

#include "execution.h"
#include "execution_policy.h"

#include "ownership.h"
#include <sycl/sycl.hpp>

namespace SPH
{
namespace execution
{
class ExecutionInstance
{
  public:
    ExecutionInstance(ExecutionInstance const &) = delete;
    void operator=(ExecutionInstance const &) = delete;

    static ExecutionInstance &getInstance()
    {
        static ExecutionInstance instance;
        return instance;
    }

    sycl::queue &getQueue()
    {
        if (!sycl_queue_)
            sycl_queue_ = makeUnique<sycl::queue>(sycl::default_selector_v);
        return *sycl_queue_;
    }

    auto getWorkGroupSize() const
    {
        return work_group_size_;
    }

    void setWorkGroupSize(size_t work_group_size)
    {
        work_group_size_ = work_group_size;
    }

    static inline sycl::nd_range<1> getUniformNdRange(size_t global_size, size_t local_size)
    {
        return {global_size % local_size ? (global_size / local_size + 1) * local_size : global_size, local_size};
    }

    inline sycl::nd_range<1> getUniformNdRange(size_t global_size) const
    {
        // sycl::nd_range is trivially-copyable, no std::move required
        return getUniformNdRange(global_size, work_group_size_);
    }

  private:
    ExecutionInstance() : work_group_size_(128), sycl_queue_() {}

    size_t work_group_size_;
    UniquePtr<sycl::queue> sycl_queue_;

} static &execution_instance = ExecutionInstance::getInstance();

} // namespace execution

/* SYCL memory transfer utilities */
template <class T>
inline T *allocateDeviceOnly(std::size_t size)
{
    return sycl::malloc_device<T>(size, execution::execution_instance.getQueue());
}

template <class T>
inline T *allocateDeviceShared(std::size_t size)
{
    return sycl::malloc_shared<T>(size, execution::execution_instance.getQueue());
}

template <class T>
inline void freeDeviceData(T *device_mem)
{
    sycl::free(device_mem, execution::execution_instance.getQueue());
}

template <class T>
inline void copyToDevice(const T *host, T *device, std::size_t size)
{
    execution::execution_instance.getQueue().memcpy(device, host, size * sizeof(T)).wait_and_throw();
}

template <class T>
inline void copyToDevice(const T &value, T *device, std::size_t size)
{
    execution::execution_instance.getQueue().fill(device, value, size).wait_and_throw();
}

template <class T>
inline void copyFromDevice(T *host, const T *device, std::size_t size)
{
    execution::execution_instance.getQueue().memcpy(host, device, size * sizeof(T)).wait_and_throw();
}

namespace execution
{
template <class LocalDynamicsType>
class Implementation<LocalDynamicsType, ParallelDevicePolicy> : public Implementation<Base>
{
    using ComputingKernel = typename LocalDynamicsType::
        template ComputingKernel<ParallelDevicePolicy>;
    UniquePtrKeeper<ComputingKernel> kernel_ptr_keeper_;

  public:
    explicit Implementation(LocalDynamicsType &local_dynamics)
        : Implementation<Base>(),
          local_dynamics_(local_dynamics), computing_kernel_(nullptr) {}
    ~Implementation()
    {
        freeDeviceData(computing_kernel_);
    }

    template <typename... Args>
    ComputingKernel *getComputingKernel(Args &&... args)
    {
        if (computing_kernel_ == nullptr)
        {
            local_dynamics_.registerComputingKernel(this, std::forward<Args>(args)...);
            computing_kernel_ = allocateDeviceOnly<ComputingKernel>(1);
            ComputingKernel *host_kernel =
                kernel_ptr_keeper_.template createPtr<ComputingKernel>(
                    ParallelDevicePolicy{}, local_dynamics_, std::forward<Args>(args)...);
            copyToDevice(host_kernel, computing_kernel_, 1);
            setUpdated();
        }

        if (!isUpdated())
        {
            overwriteComputingKernel(std::forward<Args>(args)...);
        }

        return computing_kernel_;
    }

    template <typename... Args>
    void overwriteComputingKernel(Args &&... args)
    {
        ComputingKernel *host_kernel =
            kernel_ptr_keeper_.template createPtr<ComputingKernel>(
                ParallelDevicePolicy{}, local_dynamics_, std::forward<Args>(args)...);
        copyToDevice(host_kernel, computing_kernel_, 1);
        setUpdated();
    }

  private:
    LocalDynamicsType &local_dynamics_;
    ComputingKernel *computing_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // EXECUTION_SYCL_H
