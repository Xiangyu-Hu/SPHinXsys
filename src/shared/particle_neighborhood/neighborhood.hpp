#ifndef NEIGHBORHOOD_HPP
#define NEIGHBORHOOD_HPP

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
NeighborList::ComputingKernel<ExecutionPolicy>::
    ComputingKernel(const ExecutionPolicy &ex_policy, NeighborList &neighbor_list)
    : real_particle_bound_plus_one_(neighbor_list.real_particle_bound_plus_one_),
      neighbor_index_(neighbor_list.dv_neighbor_index_->DelegatedDataField(ex_policy)),
      particle_offset_(neighbor_list.dv_particle_offset_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
template <typename FunctionOnEach>
void NeighborList::ComputingKernel<ExecutionPolicy>::
    forEachNeighbor(UnsignedInt index_i, const FunctionOnEach &function) const
{
    for (UnsignedInt n = particle_offset_[index_i]; n < particle_offset_[index_i + 1]; ++n)
    {
        function(particle_index_[n]);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBORHOOD_HPP
