#ifndef NEIGHBORHOOD_HPP
#define NEIGHBORHOOD_HPP

#include "neighborhood.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
NeighborList<ExecutionPolicy>::
    NeighborList(const ExecutionPolicy &ex_policy,
                 DiscreteVariable<UnsignedInt> *dv_neighbor_index,
                 DiscreteVariable<UnsignedInt> *dv_particle_offset)
    : neighbor_index_(dv_neighbor_index->DelegatedDataField(ex_policy)),
      particle_offset_(dv_particle_offset->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBORHOOD_HPP
