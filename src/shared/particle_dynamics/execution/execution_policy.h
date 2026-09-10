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

template <typename PolicyType>
class DecomposedExecution : public PolicyType
{
};

using MultiDevicePolicy = DecomposedExecution<SYCLDevicePolicy>;
using MultiHostPolicy = DecomposedExecution<ParallelPolicy>;
using SequencedMultiHostPolicy = DecomposedExecution<SequencedPolicy>;

inline constexpr auto seq = SequencedPolicy{};
inline constexpr auto unseq = UnsequencedPolicy{};
inline constexpr auto par_host = ParallelPolicy{};
inline constexpr auto par_unseq = ParallelUnsequencedPolicy{};
inline constexpr auto multi_device = MultiDevicePolicy{};
inline constexpr auto multi_host = MultiHostPolicy{};
inline constexpr auto seq_multi_host = SequencedMultiHostPolicy{};

#if SPHINXSYS_USE_SYCL
#if SPHINXSYS_MULTI_DEVICE
using MainExecutionPolicy = MultiDevicePolicy;
inline constexpr auto par_ck = MultiDevicePolicy{};
#else
using MainExecutionPolicy = SYCLDevicePolicy;
inline constexpr auto par_ck = SYCLDevicePolicy{};
#endif // SPHINXSYS_MULTI_DEVICE
#else
#if SPHINXSYS_MULTI_SUBDOMAIN_HOST
using MainExecutionPolicy = MultiHostPolicy;
inline constexpr auto par_ck = MultiHostPolicy{};
#else
using MainExecutionPolicy = ParallelPolicy;
inline constexpr auto par_ck = ParallelPolicy{};
#endif // SPHINXSYS_MULTI_SUBDOMAIN_HOST
#endif // SPHINXSYS_USE_SYCL

} // namespace execution
} // namespace SPH
#endif // EXECUTION_POLICY_H
