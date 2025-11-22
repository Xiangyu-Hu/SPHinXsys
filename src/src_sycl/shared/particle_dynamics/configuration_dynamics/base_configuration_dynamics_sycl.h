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
 * @file    base_configuration_dynamics_sycl.h
 * @brief   Identify data types according to execution policy.
 * @author	Xiangyu Hu
 */

#ifndef BASE_CONFIGURATION_DYNAMICS_SYCL_H
#define BASE_CONFIGURATION_DYNAMICS_SYCL_H

#include "algorithm_primitive_sycl.h"
#include "base_configuration_dynamics.h"
namespace SPH
{

class RadixSort;
template <>
struct SortMethod<ParallelDevicePolicy>
{
    typedef RadixSort type;
};

template <>
struct PlusUnsignedInt<ParallelDevicePolicy>
{
    typedef sycl::plus<UnsignedInt> type;
};
} // namespace SPH
#endif // BASE_CONFIGURATION_DYNAMICS_SYCL_H
