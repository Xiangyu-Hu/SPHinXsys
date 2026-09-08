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
#include "device_environment_sycl.h"
#include "execution_policy.h"
#include "implementation.h"
#include "ownership.h"
#include <sycl/sycl.hpp>

namespace SPH
{
namespace execution
{
/**
 * @class ExecutionInstance
 * @brief Thin facade over DeviceEnvironment resolving resources for the current device.
 * @details Everything it exposes is keyed on execution::currentSubdomainID(), which is
 *          zero unless the calling thread entered a SubdomainScope. Single device runs
 *          therefore behave exactly as before, while a multi-device fan-out gets the
 *          queue and the work group size of the device it drives, without any of the
 *          call sites having to be aware of the device at all.
 */
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
        return device_environment.getCurrentQueue();
    }

    sycl::queue &getQueue(int device_id)
    {
        return device_environment.getQueue(device_id);
    }

    auto getWorkGroupSize()
    {
        return device_environment.getWorkGroupSize(currentSubdomainID());
    }

    void setWorkGroupSize(size_t work_group_size)
    {
        device_environment.setWorkGroupSize(currentSubdomainID(), work_group_size);
    }

    static inline sycl::nd_range<1> getUniformNdRange(size_t global_size, size_t local_size)
    {
        return {global_size % local_size ? (global_size / local_size + 1) * local_size : global_size, local_size};
    }

    inline sycl::nd_range<1> getUniformNdRange(size_t global_size)
    {
        // sycl::nd_range is trivially-copyable, no std::move required
        return getUniformNdRange(global_size, getWorkGroupSize());
    }

  private:
    ExecutionInstance() = default;

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

/* Multi-device memory utilities.
 * All devices share one sycl::context, hence a pointer allocated for one device is a
 * legal argument to a copy enqueued on another device's queue. Peer-to-peer traffic
 * is therefore expressed as an ordinary memcpy on the destination queue; when the
 * link is missing the runtime stages through the host on its own. */
template <class T>
inline T *allocateDeviceOnlyOn(int device_id, std::size_t size)
{
    auto &environment = execution::device_environment;
    return sycl::malloc_device<T>(size, environment.getDevice(device_id), environment.getContext());
}

template <class T>
inline T *allocateDeviceSharedOn(int device_id, std::size_t size)
{
    auto &environment = execution::device_environment;
    return sycl::malloc_shared<T>(size, environment.getDevice(device_id), environment.getContext());
}

template <class T>
inline void freeDeviceDataOn(int device_id, T *device_mem)
{
    sycl::free(device_mem, execution::device_environment.getContext());
}

/* Copies between two devices of the shared context go through
 * execution::copyBetweenSubdomains(), declared in device_environment_sycl.h, so that
 * the exchange code stays backend independent. */

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
