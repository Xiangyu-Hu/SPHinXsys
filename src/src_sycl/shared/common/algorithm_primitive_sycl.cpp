#include "algorithm_primitive_sycl.h"

#include "sphinxsys_variable_sycl.hpp"

#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>

namespace SPH
{
//=================================================================================================//
void RadixSort::sort(const ParallelDevicePolicy &ex_policy, UnsignedInt size, UnsignedInt start_index)
{
    UnsignedInt *index_permutation = dv_index_permutation_->DelegatedData(ex_policy);
    UnsignedInt *begin = dv_sequence_->DelegatedData(ex_policy) + start_index;
    oneapi::dpl::sort_by_key(oneapi::dpl::execution::make_device_policy(execution_instance.getQueue()),
                             begin, begin + size, index_permutation);
}
//=================================================================================================//
} // namespace SPH