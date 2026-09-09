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
#include "sphinxsys_tbb.h"
#include "tbb/parallel_reduce.h"

namespace SPH
{
using namespace execution;

//----------------------------------------------------------------------
// The loop bodies are written once, on a generic loop range, and the policy
// overloads below only select between them. This is what lets the host side
// decomposition policies (MultiHostExecution<...>) reuse the same loops: a
// LoopRangeCK<ParallelMultiHostPolicy, ...> is a distinct type from a
// LoopRangeCK<ParallelPolicy, ...>, so it needs its own overload, but the loop
// itself is identical.
//----------------------------------------------------------------------
template <class LoopRangeType, class UnaryFunc>
void sequenced_particle_for(const LoopRangeType &loop_range, const UnaryFunc &unary_func)
{
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
        loop_range.computeUnit(unary_func, i);
};

template <class LoopRangeType, class UnaryFunc>
void parallel_particle_for(const LoopRangeType &loop_range, const UnaryFunc &unary_func)
{
    tbb::parallel_for(
        IndexRange(0, loop_range.LoopBound()),
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i < r.end(); ++i)
            {
                loop_range.computeUnit(unary_func, i);
            }
        },
        ap);
};

template <typename Operation, class LoopRangeType, class ReturnType, class UnaryFunc>
ReturnType sequenced_particle_reduce(const LoopRangeType &loop_range,
                                     ReturnType temp, const UnaryFunc &unary_func)
{
    Operation operation;
    ReturnType temp0 = temp;
    for (size_t i = 0; i < loop_range.LoopBound(); ++i)
    {
        temp0 = operation(temp0, loop_range.computeUnit(temp, operation, unary_func, i));
    }
    return temp0;
}

template <typename Operation, class LoopRangeType, class ReturnType, class UnaryFunc>
ReturnType parallel_particle_reduce(const LoopRangeType &loop_range,
                                    ReturnType temp, const UnaryFunc &unary_func)
{
    Operation operation;
    return tbb::parallel_reduce(
        IndexRange(0, loop_range.LoopBound()), temp,
        [&](const IndexRange &r, ReturnType temp0) -> ReturnType
        {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                temp0 = operation(temp0, loop_range.computeUnit(temp, operation, unary_func, i));
            }
            return temp0;
        },
        [&](const ReturnType &x, const ReturnType &y) -> ReturnType
        {
            return operation(x, y);
        });
};

template <class Identifier, class UnaryFunc>
void particle_for(const LoopRangeCK<SequencedPolicy, Identifier> &loop_range,
                  const UnaryFunc &unary_func)
{
    sequenced_particle_for(loop_range, unary_func);
};

template <class Identifier, class UnaryFunc>
void particle_for(const LoopRangeCK<ParallelPolicy, Identifier> &loop_range,
                  const UnaryFunc &unary_func)
{
    parallel_particle_for(loop_range, unary_func);
};

template <class Identifier, class UnaryFunc>
void particle_for(const LoopRangeCK<SequencedMultiHostPolicy, Identifier> &loop_range,
                  const UnaryFunc &unary_func)
{
    sequenced_particle_for(loop_range, unary_func);
};

template <class Identifier, class UnaryFunc>
void particle_for(const LoopRangeCK<ParallelMultiHostPolicy, Identifier> &loop_range,
                  const UnaryFunc &unary_func)
{
    parallel_particle_for(loop_range, unary_func);
};

template <typename Operation, class Identifier, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const LoopRangeCK<SequencedPolicy, Identifier> &loop_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    return sequenced_particle_reduce<Operation>(loop_range, temp, unary_func);
}

template <typename Operation, class Identifier, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const LoopRangeCK<ParallelPolicy, Identifier> &loop_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    return parallel_particle_reduce<Operation>(loop_range, temp, unary_func);
};

template <typename Operation, class Identifier, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const LoopRangeCK<SequencedMultiHostPolicy, Identifier> &loop_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    return sequenced_particle_reduce<Operation>(loop_range, temp, unary_func);
}

template <typename Operation, class Identifier, class ReturnType, class UnaryFunc>
ReturnType particle_reduce(const LoopRangeCK<ParallelMultiHostPolicy, Identifier> &loop_range,
                           ReturnType temp, const UnaryFunc &unary_func)
{
    return parallel_particle_reduce<Operation>(loop_range, temp, unary_func);
};
} // namespace SPH
#endif // PARTICLE_ITERATORS_CK_H
