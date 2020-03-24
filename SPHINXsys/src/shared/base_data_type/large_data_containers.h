#ifndef SPHINXSYS_BASE_CONTAINER_H
#define SPHINXSYS_BASE_CONTAINER_H

#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
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
	using LargeSet = tbb::concurrent_unordered_set<T>;

	typedef std::pair<std::string, std::string> StrPair;

	typedef std::pair<int, int> IndexPair;

	template<class DataType>
	using MeshData2 = DataType**;

	template<class DataType>
	using MeshData3 = DataType***;
}

#endif // SPHINXSYS_BASE_CONTAINER_H