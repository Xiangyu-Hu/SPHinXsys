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
#include "execution_policy.h"
#include "execution_queue.hpp"
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

template<class ValueType>
class DeviceRadixSort
{
    using SortablePair = std::pair<DeviceInt, ValueType>;
  public:
    /** Get the number of bits corresponding to the d-th digit of key,
     *  with each digit composed of a number of bits equal to radix_bits */
    static inline DeviceInt get_digit(DeviceInt key, DeviceInt d, DeviceInt radix_bits);

    /** Get the b-th bit of key */
    static inline DeviceInt get_bit(DeviceInt key, DeviceInt b);

    /** Group operation to compute rank, i.e. sorted position, of each work-item based on one bit.
     *  All work-items with bit = 0 will be on the first half of the ranking, while work-items with
     *  bit = 1 will be placed on the second half. */
    static SYCL_EXTERNAL DeviceInt split_count(bool bit, sycl::nd_item<1> &item);

    DeviceInt find_max_element(const DeviceInt *data, DeviceInt size, DeviceInt identity);

    void resize(DeviceInt data_size, DeviceInt radix_bits, DeviceInt workgroup_size);

    sycl::event sort_by_key(
        DeviceInt *keys, ValueType *data, DeviceInt data_size, sycl::queue &queue,
        DeviceInt workgroup_size = 256, DeviceInt radix_bits = 4);

  private:
    bool uniform_case_masking_;
    DeviceInt data_size_ = 0, radix_bits_, workgroup_size_,
           uniform_global_size_, workgroups_, radix_;
    sycl::nd_range<1> kernel_range_{0,0};
    std::unique_ptr<sycl::buffer<DeviceInt, 2>> global_buckets_, global_buckets_offsets_;
    std::unique_ptr<sycl::buffer<DeviceInt, 1>> local_buckets_offsets_buffer_;
    std::unique_ptr<sycl::buffer<SortablePair>> data_swap_buffer_, uniform_extra_swap_buffer_;
};

class BaseParticles;

struct SwapParticleDataValue
{
    template <typename DataType>
    void operator()(DataContainerAddressKeeper<StdLargeVec<DataType>> &data_keeper, size_t index_a, size_t index_b) const
    {
        for (size_t i = 0; i != data_keeper.size(); ++i)
        {
            StdLargeVec<DataType> &variable = *data_keeper[i];
            std::swap(variable[index_a], variable[index_b]);
        }
    };
};

template <typename VariableType>
struct swapParticleDeviceDataValue
{
    swapParticleDeviceDataValue() = default;

    /** Initialize swapping for device variables */
    void init(DeviceVariables &device_variables)
    {
        size_variables_ = std::get<DataTypeIndex<VariableType>::value>(device_variables).size();
        size_single_variable_ = size_variables_ ? std::get<DataTypeIndex<VariableType>::value>(device_variables)[0]->getSize()
                                                : 0;

        StdVec<DeviceVariable<VariableType> *> variables_ptr = std::get<DataTypeIndex<VariableType>::value>(device_variables);
        for (DeviceInt i = 0; i != variables_ptr.size(); ++i)
            variables_.push_back(variables_ptr[i]->VariableAddress());

        tmp_variable_ = allocateDeviceData<VariableType>(size_single_variable_);
    }

    /** Initialize swapping for extra variable not registered in DeviceVariables */
    void init(VariableType *single_variable, DeviceInt size_variable)
    {
        size_variables_ = 1;
        size_single_variable_ = size_variable;
        variables_.push_back(single_variable);
        tmp_variable_ = allocateDeviceData<VariableType>(size_variable);
    }

    ~swapParticleDeviceDataValue() { freeDeviceData(tmp_variable_); }

    /** Reorder initialized variables based on a permutation indicating where each element should be placed */
    sycl::event operator()(const DeviceInt *index_permutation, DeviceInt size) const
    {
        if (size_variables_)
            assert(size == size_single_variable_ && "Provided index permutation has different size than device variables");

        sycl::event sort_event{};
        for (DeviceInt var = 0; var < size_variables_; ++var)
        {
            auto *sortable_device_variable = variables_[var];
            auto tmp_copy_event = execution::executionQueue.getQueue().copy(sortable_device_variable, tmp_variable_, size, sort_event);
            sort_event = execution::executionQueue.getQueue().parallel_for(execution::executionQueue.getUniformNdRange(size), tmp_copy_event,
                                                                           [=, tmp_variable = tmp_variable_](sycl::nd_item<1> item)
                                                                           {
                                                                               DeviceInt i = item.get_global_id();
                                                                               if (i < size)
                                                                                   sortable_device_variable[i] = tmp_variable[index_permutation[i]];
                                                                           });
        }
        return sort_event;
    }

  private:
    DeviceInt size_variables_,                 // number of device variables registered,
        size_single_variable_;              // size of each variable
    std::vector<VariableType *> variables_; // variables to be sorted
    VariableType *tmp_variable_ = nullptr;            // temporary memory to execute sorting in parallel
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
    OperationOnDataAssemble<ParticleData, SwapParticleDataValue> swap_particle_data_value_;

  public:
    explicit SwapSortableParticleData(BaseParticles &base_particles);
    ~SwapSortableParticleData(){};

    /** the operator overload for swapping particle data.
     *  the arguments are the same with std::iter_swap
     */
    void operator()(size_t *a, size_t *b);
};

class MoveSortableParticleDeviceData
{
  protected:
    BaseParticles &base_particles_;
    bool initialized_swap_variables_;
    DeviceDataAssembleOperation<swapParticleDeviceDataValue> swap_particle_device_data_value_;
    swapParticleDeviceDataValue<DeviceInt> swap_unsorted_id_;

  public:
    explicit MoveSortableParticleDeviceData(BaseParticles &base_particles);
    /** Reorder sortable device variables based on sorting permutation*/
    void operator()(DeviceInt *index_permutation, DeviceInt size);
};

/**
 * @class ParticleSorting
 * @brief The class for sorting particle according a given sequence.
 */
class ParticleSorting
{
  protected:
    BaseParticles &base_particles_;
    DeviceInt *index_sorting_device_variables_;

    /** using pointer because it is constructed after particles. */
    SwapSortableParticleData swap_sortable_particle_data_;
    MoveSortableParticleDeviceData move_sortable_particle_device_data_;
    DeviceRadixSort<DeviceInt> device_radix_sorting;
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
    template <class ExecutionPolicy>
    inline void sortingParticleData(size_t *begin, size_t size, ExecutionPolicy execution_policy)
    {
        this->sortingParticleData(begin, size);
    }

    void sortingParticleData(DeviceInt *begin, DeviceInt size, execution::ParallelSYCLDevicePolicy execution_policy);

    /** update the reference of sorted data from unsorted data */
    virtual void updateSortedId();
    void updateSortedDeviceId() const;
};
} // namespace SPH
#endif // PARTICLE_SORTING_H