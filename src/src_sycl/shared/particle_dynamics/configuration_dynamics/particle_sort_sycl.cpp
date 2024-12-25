#include "particle_sort_sycl.h"

#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>

namespace SPH
{
//=================================================================================================//
void RadixSort::sort(const ParallelDevicePolicy &ex_policy, BaseParticles *particles)
{
    UnsignedInt *index_permutation = dv_index_permutation_->DelegatedData(ex_policy);
    UnsignedInt *sequence = dv_sequence_->DelegatedData(ex_policy);
    UnsignedInt total_real_particles = particles->TotalRealParticles();
    oneapi::dpl::sort_by_key(oneapi::dpl::execution::make_device_policy(execution_instance.getQueue()),
                             sequence, sequence + total_real_particles, index_permutation);
}
//=================================================================================================//
} // namespace SPH