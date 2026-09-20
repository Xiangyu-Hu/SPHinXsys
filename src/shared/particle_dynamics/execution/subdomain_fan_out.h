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
 * @file 	subdomain_fan_out.h
 * @brief 	Policy generic fan-out of an algorithm body over the subdomains of a run.
 * @details An algorithm's exec() wraps its body into fanOutOverSubdomains(), or into
 *          reduceOverSubdomains() when it produces a value. Under a non-decomposed
 *          policy the body simply runs once on the calling thread, so a single
 *          implementation of exec() serves every policy. Under a decomposed policy the
 *          body runs once per subdomain with that subdomain bound, which is what makes
 *          the computing kernels and the loop ranges resolve to its replica.
 *
 *          Also declares copyBetweenSubdomains(), the one operation whose meaning
 *          genuinely differs between the host and the device path: a memcpy between
 *          host allocations on one side, a queue copy between USM allocations of a
 *          shared context on the other. Isolating it here is what lets the exchange
 *          code itself be written once.
 * @author	Niki Loppi
 */

#ifndef SUBDOMAIN_FAN_OUT_H
#define SUBDOMAIN_FAN_OUT_H

#include "execution_policy.h"
#include "subdomain_runner.h"

#include <algorithm>
#include <array>
#include <type_traits>

namespace SPH
{
namespace execution
{
//----------------------------------------------------------------------
// Non-decomposed policies: the body runs once, here.
//----------------------------------------------------------------------
template <class ExecutionPolicy, class Body>
inline void fanOutOverSubdomains(const ExecutionPolicy &ex_policy, const Body &body)
{
    body();
}

template <typename Operation, class ExecutionPolicy, class ReturnType, class Body>
inline ReturnType reduceOverSubdomains(const ExecutionPolicy &ex_policy,
                                       ReturnType identity, const Body &body)
{
    return body();
}

//----------------------------------------------------------------------
// Host side decomposition, driven by the SubdomainRunner.
// The multi-device overloads live in device_environment_sycl.h.
//----------------------------------------------------------------------
/** Host decomposed policies only: the multi-device overload in device_environment_sycl.h
 *  is a non-template on MultiDevicePolicy, but the constraint keeps this one from being
 *  picked when that header is not visible at the point of instantiation. */
template <class PolicyType, class Body,
          typename = std::enable_if_t<!std::is_base_of_v<SYCLDevicePolicy, PolicyType>>>
inline void fanOutOverSubdomains(const DecomposedExecution<PolicyType> &ex_policy, const Body &body)
{
    subdomain_runner.forEachSubdomain([&](int subdomain_id)
                                      { body(); });
}

/**
 * Each subdomain reduces its own particles, and the partial results are combined in
 * subdomain order so that the result is reproducible for a given decomposition. Note
 * that it still differs from the non-decomposed result by floating point association;
 * that difference is expected and is one of the things the CPU path lets you measure
 * cheaply, ahead of interpreting a multi-GPU regression failure.
 */
template <typename Operation, class PolicyType, class ReturnType, class Body,
          typename = std::enable_if_t<!std::is_base_of_v<SYCLDevicePolicy, PolicyType>>>
inline ReturnType reduceOverSubdomains(const DecomposedExecution<PolicyType> &ex_policy,
                                       ReturnType identity, const Body &body)
{
    std::array<ReturnType, MaxSubdomains> partial_results;
    partial_results.fill(identity);
    subdomain_runner.forEachSubdomain([&](int subdomain_id)
                                      { partial_results[subdomain_id] = body(); });

    Operation operation;
    ReturnType result = identity;
    for (int subdomain_id = 0; subdomain_id < subdomain_runner.NumberOfSubdomains(); ++subdomain_id)
    {
        result = operation(result, partial_results[subdomain_id]);
    }
    return result;
}

//----------------------------------------------------------------------
// Copy between the replicas of two subdomains.
//----------------------------------------------------------------------
/** Host replicas are ordinary allocations, so the copy is a plain one and the
 *  destination subdomain argument is only there to match the device signature.
 *  std::copy rather than memcpy, since the particle data types include Eigen matrices
 *  and the assemble is instantiated for all of them. */
template <class ExecutionPolicy, class DataType>
inline void copyBetweenSubdomains(const ExecutionPolicy &ex_policy, int destination_subdomain,
                                  DataType *destination, const DataType *source, std::size_t size)
{
    std::copy(source, source + size, destination);
}

/** A decomposed policy copies exactly like its base policy: the host policies fall back to
 *  the overload above, MultiDevicePolicy reaches the queue copy of the SYCL backend. Without
 *  this forwarder the generic template would be an exact match for MultiDevicePolicy and win
 *  over the SYCLDevicePolicy overload, which needs a derived-to-base conversion. */
template <class PolicyType, class DataType>
inline void copyBetweenSubdomains(const DecomposedExecution<PolicyType> &ex_policy, int destination_subdomain,
                                  DataType *destination, const DataType *source, std::size_t size)
{
    copyBetweenSubdomains(PolicyType{}, destination_subdomain, destination, source, size);
}
} // namespace execution
} // namespace SPH
#endif // SUBDOMAIN_FAN_OUT_H
