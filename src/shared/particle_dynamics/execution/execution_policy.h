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
 * @file 	execution_policy.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Xiangyu Hu and  Fabien Pean
 */

#ifndef EXECUTION_POLICY_H
#define EXECUTION_POLICY_H

#include "subdomain_scope.h"

namespace SPH
{
namespace execution
{
class SequencedPolicy
{
};

class UnsequencedPolicy
{
};

class ParallelPolicy
{
};

class ParallelUnsequencedPolicy
{
};

class SYCLDevicePolicy
{
};

/** Tag identifying policies which fan out over the subdomains of a decomposed run. */
class MultiSubdomainTag
{
};

/**
 * @class MultiDeviceExecution
 * @brief Execution over all devices of the node, one host thread per device.
 * @details It derives from DeviceExecution so that every existing overload keyed on
 *          DeviceExecution<PolicyType> -- notably DiscreteVariable::DelegatedData and
 *          the computing kernel allocation in implementation.h -- keeps applying. The
 *          replication into per-device resources happens inside those overloads by way
 *          of execution::currentSubdomainID(), so the policy itself stays a plain tag.
 */
template <typename PolicyType>
class MultiDeviceExecution
    : public DeviceExecution<PolicyType>, public MultiSubdomainTag
{
};

using ParallelMultiDevicePolicy = MultiDeviceExecution<ParallelPolicy>;

/**
 * @class MultiHostExecution
 * @brief Domain decomposed execution entirely on the host, one subdomain at a time or
 *        one host thread per subdomain.
 * @details Its purpose is to be a debugging vehicle for MultiDeviceExecution: the two
 *          share the decomposition, the halo exchange, the migration and every fan-out
 *          point, and differ only in where the replicas live and how a copy between
 *          them is performed. A decomposition bug therefore reproduces here, under a
 *          debugger, without a GPU and without a device toolchain.
 *
 *          It derives from PolicyType so that the existing host overloads of
 *          particle_for and exclusive_scan keep applying. It does NOT derive from
 *          DeviceExecution, so DelegatedData() resolves to the host replica overload
 *          rather than to device memory.
 */
template <typename PolicyType>
class MultiHostExecution
    : public PolicyType, public MultiSubdomainTag
{
};

using ParallelMultiHostPolicy = MultiHostExecution<ParallelPolicy>;
/** Fully deterministic: one subdomain at a time, one particle at a time. */
using SequencedMultiHostPolicy = MultiHostExecution<SequencedPolicy>;

inline constexpr auto seq = SequencedPolicy{};
inline constexpr auto unseq = UnsequencedPolicy{};
inline constexpr auto par_host = ParallelPolicy{};
inline constexpr auto par_unseq = ParallelUnsequencedPolicy{};
inline constexpr auto par_device = ParallelDevicePolicy{};
inline constexpr auto seq_device = SequencedDevicePolicy{};
inline constexpr auto par_multi_device = ParallelMultiDevicePolicy{};
inline constexpr auto par_multi_host = ParallelMultiHostPolicy{};
inline constexpr auto seq_multi_host = SequencedMultiHostPolicy{};

#if SPHINXSYS_USE_SYCL
#if SPHINXSYS_MULTI_DEVICE
using MainExecutionPolicy = ParallelMultiDevicePolicy;
inline constexpr auto par_ck = ParallelMultiDevicePolicy{};
#else
using MainExecutionPolicy = ParallelDevicePolicy;
inline constexpr auto par_ck = ParallelDevicePolicy{};
#endif // SPHINXSYS_MULTI_DEVICE
#else
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
using MainExecutionPolicy = ParallelMultiHostPolicy;
inline constexpr auto par_ck = ParallelMultiHostPolicy{};
#else
using MainExecutionPolicy = ParallelPolicy;
inline constexpr auto par_ck = ParallelPolicy{};
#endif // SPHINXSYS_MULTI_SUBDOMAIN_HOST
#endif // SPHINXSYS_USE_SYCL

} // namespace execution
} // namespace SPH
#endif // EXECUTION_POLICY_H
