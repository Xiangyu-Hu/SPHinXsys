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
 * @file 	device_environment_sycl.h
 * @brief 	Multi-device SYCL environment: one queue per device on a shared context.
 * @details All devices of the node are placed in a single sycl::context. Two
 *          consequences make the multi-GPU implementation tractable:
 *          - a USM device pointer allocated for device A is a valid argument to a
 *            copy submitted on device B's queue, so peer-to-peer halo exchange is a
 *            plain queue.memcpy() without any explicit staging code path;
 *          - a kernel bundle is built once per context rather than per queue.
 *
 *          Device parallelism on the host is provided by a persistent worker pool,
 *          one thread per device. forEachSubdomain() runs the given callable on every
 *          worker with the corresponding SubdomainScope bound, and returns once all of
 *          them completed. Creating threads per call would dominate the runtime of
 *          the small kernels of an acoustic step, hence the pool.
 * @author	Niki Loppi
 */

#ifndef DEVICE_ENVIRONMENT_SYCL_H
#define DEVICE_ENVIRONMENT_SYCL_H

#include "execution_policy.h"
#include "ownership.h"
#include "subdomain_fan_out.h"
#include "subdomain_runner.h"
#include "subdomain_scope.h"

#include <array>
#include <functional>
#include <string>
#include <sycl/sycl.hpp>
#include <vector>

namespace SPH
{
namespace execution
{
/**
 * @class DeviceEnvironment
 * @brief Owns the devices, the shared context and the per-device queues.
 */
class DeviceEnvironment
{
  public:
    DeviceEnvironment(const DeviceEnvironment &) = delete;
    void operator=(const DeviceEnvironment &) = delete;

    static DeviceEnvironment &getInstance()
    {
        static DeviceEnvironment instance;
        return instance;
    }

    /** Select the devices to run on. Must be called before any device memory is
     *  allocated, i.e. before the first SPHBody is given particles.
     *  @param requested_devices number of devices to use, 0 meaning all available.
     *  Falls back to the SYCL default device when no GPU is exposed. */
    void initialize(int requested_devices = 0);

    /** Lazily initializes with all available devices on first use. */
    void ensureInitialized()
    {
        if (!initialized_)
        {
            initialize(0);
        }
    };

    int NumberOfDevices()
    {
        ensureInitialized();
        return number_of_devices_;
    };

    sycl::queue &getQueue(int device_id)
    {
        ensureInitialized();
        return *queues_[device_id];
    };

    sycl::queue &getCurrentQueue() { return getQueue(currentSubdomainID()); };

    sycl::device &getDevice(int device_id)
    {
        ensureInitialized();
        return devices_[device_id];
    };

    sycl::context &getContext()
    {
        ensureInitialized();
        return *context_;
    };

    std::size_t getWorkGroupSize(int device_id)
    {
        ensureInitialized();
        return work_group_sizes_[device_id];
    };

    void setWorkGroupSize(int device_id, std::size_t work_group_size)
    {
        work_group_sizes_[device_id] = work_group_size;
    };

    /** Whether device src can be read directly by device dst. When false the
     *  runtime still performs the copy, transparently staging through the host. */
    bool peerAccessEnabled(int src_device, int dst_device);

    /** Fan out over all devices, one host thread each, and join. */
    template <class Body>
    void forEachSubdomain(const Body &body)
    {
        ensureInitialized();
        if (number_of_devices_ == 1 || insideFanOut())
        { // single device, or a nested fan-out already bound to one device
            body(currentSubdomainID());
            return;
        }
        std::function<void(int)> wrapped = [&](int device_id)
        { body(device_id); };
        worker_pool_->run(wrapped); // the workers are permanently marked as in-fan-out
    };

    /** Wait for all queues to drain. Used at synchronization points such as I/O. */
    void synchronizeAllDevices();

    /** Human readable summary, printed once at startup. */
    std::string describe();

  private:
    DeviceEnvironment() = default;
    ~DeviceEnvironment() = default;

    void enablePeerAccess();

    bool initialized_ = false;
    int number_of_devices_ = 1;
    std::vector<sycl::device> devices_;
    UniquePtr<sycl::context> context_;
    std::vector<UniquePtr<sycl::queue>> queues_;
    std::vector<std::size_t> work_group_sizes_;
    /** peer_access_[i * number_of_devices_ + j] : device j can access device i. */
    std::vector<char> peer_access_;
    UniquePtr<SubdomainWorkerPool> worker_pool_;
};

inline DeviceEnvironment &device_environment = DeviceEnvironment::getInstance();

/** The device replicas live in USM of the shared context, so a copy between two of
 *  them is a queue operation rather than a host copy. Enqueued on the destination
 *  queue, matching the pull direction of the exchange protocol. */
template <class PolicyType, class DataType>
inline void copyBetweenSubdomains(const DeviceExecution<PolicyType> &ex_policy,
                                  int destination_subdomain, DataType *destination,
                                  const DataType *source, std::size_t size)
{
    device_environment.getQueue(destination_subdomain)
        .memcpy(destination, source, size * sizeof(DataType))
        .wait_and_throw();
}

/** Multi-device overrides of the policy generic fan-out declared in subdomain_fan_out.h. */
template <class PolicyType, class Body>
inline void fanOutOverSubdomains(const MultiDeviceExecution<PolicyType> &ex_policy, const Body &body)
{
    device_environment.forEachSubdomain([&](int device_id)
                                     { body(); });
}

/**
 * Each device reduces its own subdomain, and the partial results are combined on the
 * host. The combination is done in device order so that the result is reproducible
 * for a given decomposition; note that it still differs from the single-device result
 * by floating point association, exactly as an MPI reduction would.
 */
template <typename Operation, class PolicyType, class ReturnType, class Body>
inline ReturnType reduceOverSubdomains(const MultiDeviceExecution<PolicyType> &ex_policy,
                                    ReturnType identity, const Body &body)
{
    std::array<ReturnType, MaxSubdomains> partial_results;
    partial_results.fill(identity);
    device_environment.forEachSubdomain([&](int device_id)
                                     { partial_results[device_id] = body(); });

    Operation operation;
    ReturnType result = identity;
    for (int device_id = 0; device_id < device_environment.NumberOfDevices(); ++device_id)
    {
        result = operation(result, partial_results[device_id]);
    }
    return result;
}
} // namespace execution
} // namespace SPH
#endif // DEVICE_ENVIRONMENT_SYCL_H
