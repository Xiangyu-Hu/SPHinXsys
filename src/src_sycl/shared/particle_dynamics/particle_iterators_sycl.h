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
 * @file 	particle_iterators_sycl.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_SYCL_H
#define PARTICLE_ITERATORS_SYCL_H

#include "execution_sycl.h"
#include "particle_iterators.h"

namespace SPH
{
template <class ComputingKernelType, class ComputingKernelFunction>
inline void particle_for(Implementation<ComputingKernelType, ParallelDevicePolicy> &kernel_implementation,
                         const IndexRange &particles_range, const ComputingKernelFunction &kernel_function)
{
    auto &sycl_queue = execution_instance.getQueue();
    auto &sycl_buffer = kernel_implementation.getBuffer();
    const size_t particles_size = particles_range.size();
    sycl_queue.submit([&](sycl::handler &cgh)
                      {
            auto sycl_accessor = sycl_buffer.get_access(cgh, sycl::read_write);
            cgh.parallel_for(execution_instance.getUniformNdRange(particles_size), [=](sycl::nd_item<1> index) {
                                 if(index.get_global_id(0) < particles_size)
                                     kernel_function(index.get_global_id(0), sycl_accessor[0]);
                             }); })
        .wait_and_throw();
}
} // namespace SPH
#endif // PARTICLE_ITERATORS_SYCL_H
