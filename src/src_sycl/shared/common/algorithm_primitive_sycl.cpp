#include "algorithm_primitive_sycl.hpp"

#include "sphinxsys_variable_sycl.hpp"

#ifdef SPHINXSYS_USE_ONEDPL_SORTING
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#endif

namespace SPH
{
//=================================================================================================//
void RadixSort::sort(const ParallelDevicePolicy &ex_policy, UnsignedInt size, UnsignedInt start_index)
{
    UnsignedInt *index_permutation = dv_index_permutation_->DelegatedData(ex_policy);
    UnsignedInt *begin = dv_sequence_->DelegatedData(ex_policy) + start_index;

#ifndef SPHINXSYS_USE_ONEDPL_SORTING
    auto &sycl_queue = execution_instance.getQueue();
    device_radix_sorting.sort_by_key(begin, index_permutation, size, sycl_queue, sycl_queue.getWorkGroupSize(), 4).wait();
#else
    oneapi::dpl::sort_by_key(oneapi::dpl::execution::make_device_policy(execution_instance.getQueue()),
                             begin, begin + size, index_permutation);
#endif
}
//=================================================================================================//
} // namespace SPH