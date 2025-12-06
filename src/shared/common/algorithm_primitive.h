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
 * @file 	algorithm_primitive.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef ALGORITHM_PRIMITIVE_H
#define ALGORITHM_PRIMITIVE_H

#include "sphinxsys_variable.h"

#include <numeric>

/** this is a reformulation of tbb parallel_sort for particle data */
namespace tbb
{
namespace interface9
{
/** sorting particle */
template <typename RandomAccessIterator, typename Compare, typename SwapType>
class QuickSortRange
{
    inline size_t median_of_three(const RandomAccessIterator &array, size_t l, size_t m, size_t r) const
    {
        return comp_(array[l], array[m]) ? (comp_(array[m], array[r]) ? m : (comp_(array[l], array[r]) ? r : l))
                                         : (comp_(array[r], array[m]) ? m : (comp_(array[r], array[l]) ? r : l));
    }

    inline size_t PseudoMedianOfNine(const RandomAccessIterator &array, const QuickSortRange &range) const
    {
        size_t offset = range.size_ / 8u;
        return median_of_three(array,
                               median_of_three(array, 0, offset, offset * 2),
                               median_of_three(array, offset * 3, offset * 4, offset * 5),
                               median_of_three(array, offset * 6, offset * 7, range.size_ - 1));
    }

    size_t splitRange(QuickSortRange &range)
    {
        RandomAccessIterator array = range.begin_;
        RandomAccessIterator key0 = range.begin_;
        size_t m = PseudoMedianOfNine(array, range);
        if (m)
            swap_sortable_particle_data_(array, array + m);

        size_t i = 0;
        size_t j = range.size_;
        // Partition interval [i+1,j-1] with key *key0.
        for (;;)
        {
            __TBB_ASSERT(i < j, nullptr);
            // Loop must terminate since array[l]==*key0.
            do
            {
                --j;
                __TBB_ASSERT(i <= j, "bad ordering relation?");
            } while (comp_(*key0, array[j]));
            do
            {
                __TBB_ASSERT(i <= j, nullptr);
                if (i == j)
                    goto quick_sort_particle_partition;
                ++i;
            } while (comp_(array[i], *key0));
            if (i == j)
                goto quick_sort_particle_partition;
            swap_sortable_particle_data_(array + i, array + j);
        }
    quick_sort_particle_partition:
        // Put the partition key were it belongs
        swap_sortable_particle_data_(array + j, key0);
        // array[l..j) is less or equal to key.
        // array(j..r) is greater or equal to key.
        // array[j] is equal to key
        i = j + 1;
        size_t new_range_size = range.size_ - i;
        range.size_ = j;
        return new_range_size;
    }

  public:
    void operator=(const QuickSortRange &) = delete;
    QuickSortRange(const QuickSortRange &) = default;
    QuickSortRange() = default;

    static const size_t grainsize_ = 500;
    const Compare &comp_;
    SwapType &swap_sortable_particle_data_;
    size_t size_;
    RandomAccessIterator begin_;

    QuickSortRange(RandomAccessIterator begin,
                   size_t size, const Compare &compare, SwapType &swap_particle_data)
        : comp_(compare), swap_sortable_particle_data_(swap_particle_data),
          size_(size), begin_(begin) {}

    bool empty() const { return size_ == 0; }
    bool is_divisible() const { return size_ >= grainsize_; }

    QuickSortRange(QuickSortRange &range, split)
        : comp_(range.comp_), swap_sortable_particle_data_(range.swap_sortable_particle_data_), size_(splitRange(range))
          // +1 accounts for the pivot element, which is at its correct place
          // already and, therefore, is not included into sub-ranges.
          ,
          begin_(range.begin_ + range.size_ + 1)
    {
    }
};

/*
Description : QuickSort in Iterator format
Link        : https://stackoverflow.com/a/54976413/3547485
Ref			: http://www.cs.fsu.edu/~lacher/courses/COP4531/lectures/sorts/slide09.html
*/
template <typename RandomAccessIterator, typename Compare, typename SwapType>
RandomAccessIterator Partition(RandomAccessIterator first, RandomAccessIterator last, Compare &compare, SwapType &swap_particle_data)
{
    auto pivot = std::prev(last, 1);
    auto i = first;
    for (auto j = first; j != pivot; ++j)
    {
        // bool format
        if (compare(*j, *pivot))
        {
            swap_particle_data(i++, j);
        }
    }
    swap_particle_data(i, pivot);
    return i;
}

template <typename RandomAccessIterator, typename Compare, typename SwapType>
void SerialQuickSort(RandomAccessIterator first, RandomAccessIterator last, Compare &compare, SwapType &swap_particle_data)
{
    if (std::distance(first, last) > 1)
    {
        RandomAccessIterator bound = Partition(first, last, compare, swap_particle_data);
        SerialQuickSort(first, bound, compare, swap_particle_data);
        SerialQuickSort(bound + 1, last, compare, swap_particle_data);
    }
}

/*
Description : InsertionSort in Iterator format
Link        : http://www.codecodex.com/wiki/Insertion_sort
*/

template <typename RandomAccessIterator, typename Compare, typename SwapType>
void InsertionSort(RandomAccessIterator First, RandomAccessIterator Last, Compare &compare, SwapType &swap_particle_data)
{
    RandomAccessIterator min = First;
    for (RandomAccessIterator i = First + 1; i < Last; ++i)
        if (compare(*i, *min))
            min = i;

    swap_particle_data(First, min);
    while (++First < Last)
        for (RandomAccessIterator j = First; compare(*j, *(j - 1)); --j)
            swap_particle_data((j - 1), j);
}

/** Body class used to sort elements in a range that is smaller than the grainsize. */
template <typename RandomAccessIterator, typename Compare, typename SwapType>
struct QuickSortBody
{
    void operator()(const QuickSortRange<RandomAccessIterator, Compare, SwapType> &range) const
    {
        SerialQuickSort(range.begin_, range.begin_ + range.size_, range.comp_, range.swap_sortable_particle_data_);
    }
};
} // namespace interface9
} // namespace tbb

namespace SPH
{
struct CompareSequence
{
    bool operator()(const size_t &x, const size_t &y) const
    {
        return x < y;
    };
};

class QuickSort
{
    class SwapIndex
    {
        UnsignedInt *sequence_;
        UnsignedInt *index_permutation_;

      public:
        SwapIndex(UnsignedInt *sequence, UnsignedInt *index_permutation);
        ~SwapIndex() {};

        void operator()(UnsignedInt *a, UnsignedInt *b);
    };

  public:
    template <class ExecutionPolicy>
    explicit QuickSort(const ExecutionPolicy &ex_policy,
                       DiscreteVariable<UnsignedInt> *dv_sequence,
                       DiscreteVariable<UnsignedInt> *dv_index_permutation)
        : sequence_(dv_sequence->DelegatedData(ex_policy)),
          index_permutation_(dv_index_permutation->DelegatedData(ex_policy)),
          swap_index_(sequence_, index_permutation_), compare_(),
          quick_sort_range_(sequence_, 0, compare_, swap_index_),
          quick_sort_body_(){};
    void sort(const ParallelPolicy &ex_policy, UnsignedInt size, UnsignedInt start_index = 0);

  protected:
    UnsignedInt *sequence_;
    UnsignedInt *index_permutation_;
    SwapIndex swap_index_;
    CompareSequence compare_;
    tbb::interface9::QuickSortRange<UnsignedInt *, CompareSequence, SwapIndex> quick_sort_range_;
    tbb::interface9::QuickSortBody<UnsignedInt *, CompareSequence, SwapIndex> quick_sort_body_;
};

template <typename T, typename Op>
T exclusive_scan(const SequencedPolicy &seq_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    UnsignedInt scan_size = d_size - 1;
    std::exclusive_scan(first, first + d_size, d_first, T{0}, op);
    return d_first[scan_size];
}

template <typename T, typename Op>
T exclusive_scan(const ParallelPolicy &par_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    // Exclusive scan is the same as inclusive, but shifted by one
    UnsignedInt scan_size = d_size - 1;
    d_first[0] = T{0};
    using range_type = tbb::blocked_range<UnsignedInt>;
    tbb::parallel_scan(
        range_type(0, scan_size), d_first[0],
        [=](const range_type &r, T sum, bool is_final_scan) -> T
        {
            T tmp = sum;
            for (UnsignedInt i = r.begin(); i < r.end(); ++i)
            {
                tmp = op(tmp, first[i]);
                if (is_final_scan)
                {
                    d_first[i + 1] = tmp;
                }
            }
            return tmp;
        },
        [&](const T &a, const T &b)
        {
            return op(a, b);
        });
    return d_first[scan_size];
}

template <class LocalDynamicsFunction>
inline void generic_for(const SequencedPolicy &seq, const IndexRange &index_range,
                         const LocalDynamicsFunction &local_dynamics_function)
{
    for (size_t i = index_range.begin(); i < index_range.end(); ++i)
        local_dynamics_function(i);
};

template <class LocalDynamicsFunction>
inline void generic_for(const ParallelPolicy &par_host, const IndexRange &particles_range,
                         const LocalDynamicsFunction &local_dynamics_function)
{
    tbb::parallel_for(
        particles_range,
        [&](const IndexRange &r)
        {
            for (size_t i = r.begin(); i < r.end(); ++i)
            {
                local_dynamics_function(i);
            }
        },
        ap);
};
} // namespace SPH

#endif // ALGORITHM_PRIMITIVE_H
