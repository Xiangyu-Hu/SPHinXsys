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
 * @file particle_sorting.h
 * @brief Here gives the classes for particle sorting.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_SORTING_H
#define PARTICLE_SORTING_H

#include "base_data_package.h"
#include "sph_data_containers.h"

/** this is a reformulation of tbb parallel_sort for particle data */
namespace tbb
{
namespace interface9
{
/** sorting particle */
template <typename RandomAccessIterator, typename Compare, typename SwapType>
class QuickSortParticleRange
{
    inline size_t median_of_three(const RandomAccessIterator &array, size_t l, size_t m, size_t r) const
    {
        return comp_(array[l], array[m]) ? (comp_(array[m], array[r]) ? m : (comp_(array[l], array[r]) ? r : l))
                                         : (comp_(array[r], array[m]) ? m : (comp_(array[r], array[l]) ? r : l));
    }

    inline size_t PseudoMedianOfNine(const RandomAccessIterator &array, const QuickSortParticleRange &range) const
    {
        size_t offset = range.size_ / 8u;
        return median_of_three(array,
                               median_of_three(array, 0, offset, offset * 2),
                               median_of_three(array, offset * 3, offset * 4, offset * 5),
                               median_of_three(array, offset * 6, offset * 7, range.size_ - 1));
    }

    size_t splitRange(QuickSortParticleRange &range)
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
    void operator=(const QuickSortParticleRange &) = delete;
    QuickSortParticleRange(const QuickSortParticleRange &) = default;
    QuickSortParticleRange() = default;

    static const size_t grainsize_ = 500;
    const Compare &comp_;
    SwapType &swap_sortable_particle_data_;
    size_t size_;
    RandomAccessIterator begin_;

    QuickSortParticleRange(RandomAccessIterator begin,
                           size_t size, const Compare &compare, SwapType &swap_particle_data)
        : comp_(compare), swap_sortable_particle_data_(swap_particle_data),
          size_(size), begin_(begin) {}

    bool empty() const { return size_ == 0; }
    bool is_divisible() const { return size_ >= grainsize_; }

    QuickSortParticleRange(QuickSortParticleRange &range, split)
        : comp_(range.comp_), swap_sortable_particle_data_(range.swap_sortable_particle_data_), size_(splitRange(range))
          // +1 accounts for the pivot element, which is at its correct place
          // already and, therefore, is not included into subranges.
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
struct QuickSortParticleBody
{
    void operator()(const QuickSortParticleRange<RandomAccessIterator, Compare, SwapType> &range) const
    {
        SerialQuickSort(range.begin_, range.begin_ + range.size_, range.comp_, range.swap_sortable_particle_data_);
    }
};
} // namespace interface9
} // namespace tbb
/**
 * SPH implementation.
 */
namespace SPH
{
class BaseParticles;

template <typename VariableType>
struct swapParticleDataValue
{
    void operator()(ParticleData &particle_data, size_t index_a, size_t index_b) const
    {
        constexpr int type_index = DataTypeIndex<VariableType>::value;

        StdVec<StdLargeVec<VariableType> *> variables = std::get<type_index>(particle_data);
        for (size_t i = 0; i != variables.size(); ++i)
        {
            StdLargeVec<VariableType> &variable = *variables[i];
            std::swap(variable[index_a], variable[index_b]);
        }
    };
};

/**
 * @class CompareParticleSequence
 * @brief compare the sequence of two particles
 */
struct CompareParticleSequence
{
    bool operator()(const size_t &x, const size_t &y) const
    {
        return x < y;
    };
};

/**
 * @class SwapSortableParticleData
 * @brief swap sortable particle data according to a sequence
 */
class SwapSortableParticleData
{
  protected:
    StdLargeVec<size_t> &sequence_;
    StdLargeVec<size_t> &unsorted_id_;
    ParticleData &sortable_data_;
    DataAssembleOperation<swapParticleDataValue> swap_particle_data_value_;

  public:
    explicit SwapSortableParticleData(BaseParticles &base_particles);
    ~SwapSortableParticleData(){};

    /** the operator overload for swapping particle data.
     *  the arguments are the same with std::iter_swap
     */
    void operator()(size_t *a, size_t *b);
};

/**
 * @class ParticleSorting
 * @brief The class for sorting particle according a given sequence.
 */
class ParticleSorting
{
  protected:
    BaseParticles &base_particles_;

    /** using pointer because it is constructed after particles. */
    SwapSortableParticleData swap_sortable_particle_data_;
    CompareParticleSequence compare_;
    tbb::interface9::QuickSortParticleRange<
        size_t *, CompareParticleSequence, SwapSortableParticleData>
        quick_sort_particle_range_;
    tbb::interface9::QuickSortParticleBody<
        size_t *, CompareParticleSequence, SwapSortableParticleData>
        quick_sort_particle_body_;

  public:
    // the construction is before particles
    explicit ParticleSorting(BaseParticles &base_particles);
    virtual ~ParticleSorting(){};
    /** sorting particle data according to the cell location of particles */
    virtual void sortingParticleData(size_t *begin, size_t size);
    /** update the reference of sorted data from unsorted data */
    virtual void updateSortedId();
};
} // namespace SPH
#endif // PARTICLE_SORTING_H