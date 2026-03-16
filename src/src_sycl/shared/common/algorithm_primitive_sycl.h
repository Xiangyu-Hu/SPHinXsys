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
 * @file algorithm_primitive_sycl.h
 * @brief Here gives the classes for particle sorting.
 * @author Xiangyu Hu
 */

#ifndef ALGORITHM_PRIMITIVE_SYCL_H
#define ALGORITHM_PRIMITIVE_SYCL_H

#include "implementation_sycl.h"
#include "sphinxsys_variable.h"

namespace SPH
{
using namespace execution;
template <class DataType>
class DeviceRadixSort
{
    using SortablePair = std::pair<UnsignedInt, DataType>;

  public:
    /** Get the number of bits corresponding to the d-th digit of key,
     *  with each digit composed of a number of bits equal to radix_bits */
    static inline UnsignedInt get_digit(UnsignedInt key, UnsignedInt d, UnsignedInt radix_bits);

    /** Get the b-th bit of key */
    static inline UnsignedInt get_bit(UnsignedInt key, UnsignedInt b);

    /** Group operation to compute rank, i.e. sorted position, of each work-item based on one bit.
     *  All work-items with bit = 0 will be on the first half of the ranking, while work-items with
     *  bit = 1 will be placed on the second half. */
    static UnsignedInt split_count(bool bit, sycl::nd_item<1> &item);

    UnsignedInt find_max_element(const UnsignedInt *data, UnsignedInt size, UnsignedInt identity);

    void resize(UnsignedInt data_size, UnsignedInt radix_bits, UnsignedInt workgroup_size);

    sycl::event sort_by_key(
        UnsignedInt *keys, DataType *data, UnsignedInt data_size, sycl::queue &queue,
        UnsignedInt workgroup_size = 256, UnsignedInt radix_bits = 4);

  private:
    bool uniform_case_masking_;
    UnsignedInt data_size_ = 0;
    UnsignedInt radix_bits_, workgroup_size_, uniform_global_size_, workgroups_, radix_;
    sycl::nd_range<1> kernel_range_{0, 0};
    std::unique_ptr<sycl::buffer<UnsignedInt, 2>> global_buckets_, global_buckets_offsets_;
    std::unique_ptr<sycl::buffer<UnsignedInt, 1>> local_buckets_offsets_buffer_;
    std::unique_ptr<sycl::buffer<SortablePair>> data_swap_buffer_, uniform_extra_swap_buffer_;
};

class RadixSort
{
  public:
    template <class ExecutionPolicy>
    explicit RadixSort(const ExecutionPolicy &ex_policy,
                       DiscreteVariable<UnsignedInt> *dv_sequence,
                       DiscreteVariable<UnsignedInt> *dv_index_permutation)
        : dv_sequence_(dv_sequence), dv_index_permutation_(dv_index_permutation){};
    void sort(const ParallelDevicePolicy &ex_policy, UnsignedInt size, UnsignedInt start_index = 0);

  protected:
    DiscreteVariable<UnsignedInt> *dv_sequence_;
    DiscreteVariable<UnsignedInt> *dv_index_permutation_;
    DeviceRadixSort<UnsignedInt> device_radix_sorting;
};

template <typename T, typename Op>
T exclusive_scan(const ParallelDevicePolicy &par_policy, T *first, T *d_first, UnsignedInt d_size, Op op)
{
    execution_instance.getQueue()
        .submit([=](sycl::handler &cgh)
                { cgh.parallel_for(
                      execution_instance.getUniformNdRange(execution_instance.getWorkGroupSize()),
                      [=](sycl::nd_item<1> item)
                      {
                          if (item.get_group_linear_id() == 0)
                          {
                              sycl::joint_exclusive_scan(
                                  item.get_group(), first, first + d_size, d_first, T{0}, op);
                          }
                      }); })
        .wait_and_throw();

    UnsignedInt scan_size = d_size - 1;
    T last_value;
    copyFromDevice(&last_value, d_first + scan_size, 1);
    return last_value;
}

template <class UnaryFunc>
void generic_for(const ParallelDevicePolicy &par_device,
                 const IndexRange &index_range,
                 const UnaryFunc &unary_func)
{
    auto &sycl_queue = execution_instance.getQueue();
    const size_t loop_bound = index_range.size();
    sycl_queue.submit([&](sycl::handler &cgh)
                      { cgh.parallel_for(execution_instance.getUniformNdRange(loop_bound), [=](sycl::nd_item<1> index)
                                         {
                                 if(index.get_global_id(0) < loop_bound)
                                     unary_func(index.get_global_id(0)); }); })
        .wait_and_throw();
}
} // namespace SPH
#endif // ALGORITHM_PRIMITIVE_SYCL_H
