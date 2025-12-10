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
 * @file 	implementation_sycl.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Alberto Guarnieri and Xiangyu Hu
 */

#ifndef IMPLEMENTATION_SYCL_H
#define IMPLEMENTATION_SYCL_H

#include "device_copyable_variable.h"
#include "execution_policy.h"
#include "implementation.h"
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
        {
            sycl_queue_ = makeUnique<sycl::queue>(sycl::default_selector_v);
            auto device = sycl_queue_->get_device();
            unsigned long max_workgroup_size = device.get_info<sycl::info::device::max_work_group_size>();
            work_group_size_ = SMIN(max_workgroup_size, 64UL);
        }
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
inline T *allocateHostStaging(std::size_t size)
{
    return sycl::malloc_host<T>(size, execution::execution_instance.getQueue());
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
template <class ComputingKernelType>
inline ComputingKernelType *allocateComputingKernelOnDevice()
{
    return allocateDeviceOnly<ComputingKernelType>(1);
}

template <class ComputingKernelType>
inline void copyComputingKernelToDevice(ComputingKernelType *host_kernel,
                                        ComputingKernelType *device_kernel)
{
    copyToDevice(host_kernel, device_kernel, 1);
}

template <class ComputingKernelType>
inline void freeComputingKernelOnDevice(ComputingKernelType *device_kernel)
{
    freeDeviceData(device_kernel);
}
} // namespace execution
} // namespace SPH
#endif // IMPLEMENTATION_SYCL_H
