#ifndef TBB_MIDDLE_H
#define TBB_MIDDLE_H
#pragma GCC system_header
#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#pragma warning(push, 0)

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
#pragma warning(pop)
#endif //TBB_MIDDLE_H