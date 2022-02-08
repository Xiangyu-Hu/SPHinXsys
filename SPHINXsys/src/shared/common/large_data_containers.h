/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
#ifndef LARGE_DATA_CONTAINER_H
#define LARGE_DATA_CONTAINER_H

#include "tbb/tbb.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/blocked_range3d.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/tick_count.h"
#include "tbb/scalable_allocator.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/cache_aligned_allocator.h"

#include <array>

using namespace tbb;
static tbb::affinity_partitioner ap;

namespace SPH {

	template <typename T>
	using LargeVec = tbb::concurrent_vector<T>;

	template <typename T>
	using StdLargeVec = std::vector<T, cache_aligned_allocator<T>>;

	template <typename T>
	using StdVec = std::vector<T>;

	template <typename T>
	using DoubleVec = std::vector<std::vector<T>>;

	template <typename T>
	using TripleVec = std::vector<std::vector<std::vector<T>>>;
}

#endif //LARGE_DATA_CONTAINER_H
