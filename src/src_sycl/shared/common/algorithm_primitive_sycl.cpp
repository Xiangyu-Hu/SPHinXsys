#include "algorithm_primitive_sycl.hpp"

#include "sphinxsys_variable_sycl.hpp"

namespace SPH
{
//=================================================================================================//
void RadixSort::sort(const ParallelDevicePolicy &ex_policy, UnsignedInt size, UnsignedInt start_index)
{
    UnsignedInt *index_permutation = dv_index_permutation_->DelegatedData(ex_policy);
    UnsignedInt *begin = dv_sequence_->DelegatedData(ex_policy) + start_index;

#if !SPHINXSYS_USE_ONEDPL
    auto &sycl_queue = execution_instance.getQueue();
    auto work_group_size = execution_instance.getWorkGroupSize();
    device_radix_sorting.sort_by_key(begin, index_permutation, size, sycl_queue, work_group_size, 4)
        .wait_and_throw();
#else
    oneapi::dpl::execution::device_policy<> policy(execution_instance.getQueue());
    oneapi::dpl::sort_by_key(policy, begin, begin + size, index_permutation);
    policy.queue().wait_and_throw();
#endif
}
//=================================================================================================//
} // namespace SPH