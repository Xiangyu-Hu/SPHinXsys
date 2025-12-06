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
    device_radix_sorting.sort_by_key(begin, index_sorting_device_variables_, size, execution::executionQueue.getQueue(), 512, 4).wait();
#else
    execution::executionQueue.getQueue().parallel_for(execution::executionQueue.getUniformNdRange(size),
                                                      [=, index_sorting = index_sorting_device_variables_](sycl::nd_item<1> it)
                                                      {
                                                          DeviceInt i = it.get_global_id(0);
                                                          if (i < size)
                                                              index_sorting[i] = i;
                                                      })
        .wait();

    oneapi::dpl::sort_by_key(oneapi::dpl::execution::make_device_policy(execution_instance.getQueue()),
                             begin, begin + size, index_permutation);
#endif
}
//=================================================================================================//
} // namespace SPH