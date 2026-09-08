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
 * @file 	particle_iterators_multi_device_sycl.h
 * @brief 	Particle iterators for the multi-device (multi-GPU) execution policy.
 * @details The loop range of a multi-device run is not one range but one per device,
 *          each bound to the subdomain replica of that device. The iterators here do
 *          not fan out themselves: they are entered from within a device thread
 *          (see fanOutOverSubdomains) and simply pick the range of the current device,
 *          then reuse the single-device iterators verbatim. Keeping the fan-out at
 *          the algorithm level, rather than at the loop level, means a device thread
 *          keeps its device binding across the whole exec(), which is what lets the
 *          computing kernels resolve to the right replica.
 * @author	Niki Loppi
 */

#ifndef PARTICLE_ITERATORS_MULTI_DEVICE_SYCL_H
#define PARTICLE_ITERATORS_MULTI_DEVICE_SYCL_H

#include "device_environment_sycl.h"
#include "particle_iterators_sycl.h"

namespace SPH
{
/**
 * @class LoopRangeCK<ParallelMultiDevicePolicy, Identifier>
 * @brief The loop range of the device bound to the calling thread.
 * @details The CK algorithms construct their loop range inside exec(), which under the
 *          multi-device policy already runs on a device thread. The range therefore
 *          only has to resolve to the current device, and it does so by deriving from
 *          the single-device range: its constructor calls DelegatedData(), which is
 *          keyed on currentSubdomainID(). Loop bound, particle list and cell offsets are
 *          consequently those of the subdomain, not of the global particle set.
 *
 *          Invariant: a loop range must be constructed inside the fan-out. Building
 *          one on the host thread and reusing it inside the device threads would bind
 *          every device to the replica of device 0.
 *
 *          Deriving also makes the single-device particle_for/particle_reduce
 *          overloads of particle_iterators_sycl.h apply unchanged. particle_reduce
 *          then yields this device's partial value, and reduceOverSubdomains() in the
 *          enclosing ReduceDynamicsCK::exec() combines the partials across devices.
 */
template <class Identifier>
class LoopRangeCK<ParallelMultiDevicePolicy, Identifier>
    : public LoopRangeCK<ParallelDevicePolicy, Identifier>
{
  public:
    template <typename... Args>
    explicit LoopRangeCK(Args &&...args)
        : LoopRangeCK<ParallelDevicePolicy, Identifier>(std::forward<Args>(args)...) {};
};
} // namespace SPH
#endif // PARTICLE_ITERATORS_MULTI_DEVICE_SYCL_H
