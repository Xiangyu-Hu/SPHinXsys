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
 * @file    base_configuration_dynamics.h
 * @brief   TBD.
 * @author	Xiangyu Hu
 */

#ifndef BASE_CONFIGURATION_DYNAMICS_H
#define BASE_CONFIGURATION_DYNAMICS_H

#include "base_data_package.h"
#include "execution_policy.h"

#include <boost/atomic/atomic_ref.hpp>
#include <functional>

namespace SPH
{
using namespace execution;

template <typename... T>
struct AtomicUnsignedIntRef;

template <>
struct AtomicUnsignedIntRef<SequencedPolicy>
{
    typedef UnsignedInt &type;
};

template <>
struct AtomicUnsignedIntRef<ParallelPolicy>
{
    typedef boost::atomic_ref<UnsignedInt> type;
};

template <typename... T>
struct PlusUnsignedInt;

template <class ExecutionPolicy>
struct PlusUnsignedInt<ExecutionPolicy>
{
    typedef std::plus<UnsignedInt> type;
};
} // namespace SPH
#endif // BASE_CONFIGURATION_DYNAMICS_H
