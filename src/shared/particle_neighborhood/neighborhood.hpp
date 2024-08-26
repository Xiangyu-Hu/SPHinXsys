#ifndef NEIGHBORHOOD_HPP
#define NEIGHBORHOOD_HPP

#include "neighborhood.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
NeighborList::NeighborList(const ExecutionPolicy &ex_policy,
                           DiscreteVariable<UnsignedInt> &dv_neighbor_index,
                           DiscreteVariable<UnsignedInt> &dv_particle_offset)
    : neighbor_index_(dv_neighbor_index.DelegatedDataField(ex_policy)),
      particle_offset_(dv_particle_offset.DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <typename FunctionOnEach>
void NeighborList::forEachNeighbor(UnsignedInt index_i, const FunctionOnEach &function) const
{
    for (UnsignedInt n = particle_offset_[index_i]; n < particle_offset_[index_i + 1]; ++n)
    {
        function(neighbor_index_[n]);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBORHOOD_HPP
