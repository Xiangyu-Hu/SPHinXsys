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
                       DiscreteVariable<UnsignedInt> *dv_index_permutation);
    void sort(const ParallelDevicePolicy &ex_policy, UnsignedInt size);

  protected:
    DiscreteVariable<UnsignedInt> *dv_sequence_;
    DiscreteVariable<UnsignedInt> *dv_index_permutation_;
};
} // namespace SPH
#endif // ALGORITHM_PRIMITIVE_SYCL_H
