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
            sycl_queue_ = makeUnique<sycl::queue>(sycl::gpu_selector_v);
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
    ExecutionInstance() : work_group_size_(32), sycl_queue_() {}

    size_t work_group_size_;
    UniquePtr<sycl::queue> sycl_queue_;

} static &execution_instance = ExecutionInstance::getInstance();

class ExecutionEvent
{
  public:
    ExecutionEvent() = default;
    ExecutionEvent(const sycl::event &event);
    ExecutionEvent(const std::vector<sycl::event> &event);

    ExecutionEvent(const ExecutionEvent &executionEvent) = default;
    ExecutionEvent(ExecutionEvent &&executionEvent) = default;

    ExecutionEvent &operator=(const ExecutionEvent &event) = default;
    ExecutionEvent &operator=(ExecutionEvent &&event) = default;

    ExecutionEvent &operator=(sycl::event event);

    const std::vector<sycl::event> &getEventList() const;

    ExecutionEvent &add(sycl::event event);
    ExecutionEvent &add(const ExecutionEvent &event);

    void wait();
    ExecutionEvent &then(std::function<void()> &&func,
                         std::optional<std::reference_wrapper<ExecutionEvent>> host_event = {});

  private:
    std::vector<sycl::event> event_list_;
};

template <class ComputingKernelType>
class Implementation<ComputingKernelType, ParallelDevicePolicy>
{
    using DeviceComputingKernelType = sycl::buffer<ComputingKernelType>;

  public:
    Implementation() = default;

    explicit Implementation(ComputingKernelType *computing_kernel)
        : delegated_kernel_(computing_kernel, 1) {}

    DeviceComputingKernelType &getBuffer()
    {
        return delegated_kernel_;
    }

  private:
    DeviceComputingKernelType delegated_kernel_;
};
} // namespace execution
} // namespace SPH
#endif // EXECUTION_SYCL_H
