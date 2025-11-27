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
