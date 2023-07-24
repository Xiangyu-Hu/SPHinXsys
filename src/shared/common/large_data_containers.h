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
 * @file 	large_data-container.h
 * @brief 	Data container for large vector, e.g. particle data.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef LARGE_DATA_CONTAINER_H
#define LARGE_DATA_CONTAINER_H

#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/blocked_range3d.h"
#include "tbb/cache_aligned_allocator.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/scalable_allocator.h"
#include "tbb/tick_count.h"

#include <array>

namespace SPH
{

static tbb::affinity_partitioner ap;
typedef tbb::blocked_range<size_t> IndexRange;
typedef tbb::blocked_range2d<size_t> IndexRange2d;
typedef tbb::blocked_range3d<size_t> IndexRange3d;

typedef tbb::tick_count TickCount;
typedef tbb::tick_count::interval_t TimeInterval;

template <typename T>
using ConcurrentVec = tbb::concurrent_vector<T>;

template <typename T>
using StdLargeVec = std::vector<T, tbb::cache_aligned_allocator<T>>;

template <typename T>
using StdVec = std::vector<T>;

template <typename T>
using BiVector = std::vector<std::vector<T>>;

template <typename T>
using TriVector = std::vector<std::vector<std::vector<T>>>;
} // namespace SPH

#endif // LARGE_DATA_CONTAINER_H
