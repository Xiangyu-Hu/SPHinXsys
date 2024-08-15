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
 * @file    update_cell_linked_list_sycl.h
 * @brief   TBD
 * @author	Xiangyu Hu
 */

#ifndef UPDATE_CELL_LINKED_LIST_SYCL_H
#define UPDATE_CELL_LINKED_LIST_SYCL_H

#include "execution_sycl.h"
#include "update_cell_linked_list.h"

namespace SPH
{

template <>
struct AtomicUnsignedIntRef<ParallelDevicePolicy>
{
    typedef sycl::atomic_ref<
        UnsignedInt, sycl::memory_order_relaxed, sycl::memory_scope_device,
        sycl::access::address_space::global_space>
        type;
};

template <typename T, typename Op>
T exclusive_scan(const ParallelDevicePolicy &par_policy, T *first, T *last, T *d_first, Op op)
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
                                  item.get_group(), first, last, d_first, T{0}, op);
                          }
                      }); })
        .wait_and_throw();

    UnsignedInt scan_size = last - first - 1;
    UnsignedInt last_value;
    copyFromDevice(&last_value, d_first + scan_size, 1);
    return last_value;
}

} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_SYCL_H
