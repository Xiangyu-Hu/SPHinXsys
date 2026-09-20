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
 * @file 	particle_iterators_ck.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_CK_H
#define PARTICLE_ITERATORS_CK_H

#include "loop_range.h"
#include "subdomain_fan_out.h"
#include "sphinxsys_tbb.h"
#include "tbb/parallel_reduce.h"

namespace SPH
{
using namespace execution;

template <class Identifier, class KernelImplementationType>
void particle_for(const LoopRangeCK<SequencedPolicy, Identifier> &loop_range,
                  KernelImplementationType &implementation, Real dt)
{
    auto kernel = implementation.getComputingKernel();
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
        loop_range.computeUnit([=](size_t i)
                               { kernel->compute(i, dt); }, i);
};

template <class Identifier, class KernelImplementationType>
void particle_for(const LoopRangeCK<ParallelPolicy, Identifier> &loop_range,
                  KernelImplementationType &implementation, Real dt)
{
    auto kernel = implementation.getComputingKernel();
    tbb::parallel_for(
        IndexRange(0, loop_range.LoopBound()),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i < r.end(); ++i)
            {
                loop_range.computeUnit(
                    [=](size_t i)
                    { kernel->compute(i, dt); }, i);
            }
        },
        ap);
};

template <typename Operation, class Identifier, class ReturnType, class KernelImplementationType>
ReturnType particle_reduce(const LoopRangeCK<SequencedPolicy, Identifier> &loop_range,
                           ReturnType temp, KernelImplementationType &implementation, Real dt)
{
    auto reduce_kernel = implementation.getComputingKernel();
    Operation operation;
    ReturnType temp0 = temp;
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
    {
        temp0 = operation(
            temp0, loop_range.computeUnit(
                       temp, operation,
                       [=](size_t i)
                       { return reduce_kernel->reduce(i, dt); }, i));
    }
    return temp0;
}

template <typename Operation, class Identifier, class ReturnType, class KernelImplementationType>
ReturnType particle_reduce(const LoopRangeCK<ParallelPolicy, Identifier> &loop_range,
                           ReturnType temp, KernelImplementationType &implementation, Real dt)
{
    auto reduce_kernel = implementation.getComputingKernel();
    Operation operation;
    return tbb::parallel_reduce(
        IndexRange(0, loop_range.LoopBound()), temp,
        [&](const IndexRange &r, ReturnType temp0) -> ReturnType
        {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = operation(temp0, loop_range.computeUnit(
                        temp, operation,                         
                        [=](size_t i)
                       { return reduce_kernel->reduce(i, dt); }, i));
				}
				return temp0; },
        [&](const ReturnType &x, const ReturnType &y) -> ReturnType
        {
            return operation(x, y);
        });
};

template <typename Operation, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const ParallelPolicy &par, const IndexRange &particles_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    Operation operation;
    return tbb::parallel_reduce(
        particles_range, temp,
        [&](const IndexRange &r, ReturnType temp0) -> ReturnType
        {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = operation(temp0, unary_func(i));
				}
				return temp0; },
        [&](const ReturnType &x, const ReturnType &y) -> ReturnType
        {
            return operation(x, y);
        });
};

//----------------------------------------------------------------------
// Decomposed policies: the fan-out over the subdomains happens here, at the loop level.
// Each subdomain builds its own range and fetches its own computing kernel inside the
// fan-out (getComputingKernel() resolves on currentSubdomainID()), then runs the loop of
// the base policy unchanged. A fan-out is a barrier over the subdomains, so consecutive
// particle_for() calls of one algorithm are ordered across subdomains, which is what
// lets a halo refresh sit between any two of them.
//----------------------------------------------------------------------
template <class PolicyType, class Identifier, class KernelImplementationType>
void particle_for(const LoopRangeCK<DecomposedExecution<PolicyType>, Identifier> &loop_range,
                  KernelImplementationType &implementation, Real dt)
{
    fanOutOverSubdomains(DecomposedExecution<PolicyType>{}, [&]()
                         { particle_for(loop_range.onCurrentSubdomain(), implementation, dt); });
};

/** Each subdomain reduces its owned particles; the partial results are combined in
 *  subdomain order by reduceOverSubdomains(), so the result is reproducible for a given
 *  decomposition. Halo particles lie past the owned range and are never counted. */
template <typename Operation, class PolicyType, class Identifier, class ReturnType, class KernelImplementationType>
ReturnType particle_reduce(const LoopRangeCK<DecomposedExecution<PolicyType>, Identifier> &loop_range,
                           ReturnType temp, KernelImplementationType &implementation, Real dt)
{
    return reduceOverSubdomains<Operation>(
        DecomposedExecution<PolicyType>{}, temp, [&]()
        { return particle_reduce<Operation>(loop_range.onCurrentSubdomain(), temp, implementation, dt); });
};
} // namespace SPH
#endif // PARTICLE_ITERATORS_CK_H
