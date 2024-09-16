#ifndef PARTICLE_SORT_SYCL_HPP
#define PARTICLE_SORT_SYCL_HPP

#include "particle_sort_sycl.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
RadixSort::RadixSort(const ExecutionPolicy &ex_policy,
                     DiscreteVariable<UnsignedInt> *dv_sequence,
                     DiscreteVariable<UnsignedInt> *dv_index_permutation)
    : dv_sequence_(dv_sequence), dv_index_permutation_(dv_index_permutation) {}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SORT_SYCL_HPP
