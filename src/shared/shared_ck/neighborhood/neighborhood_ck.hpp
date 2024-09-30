#ifndef NEIGHBORHOOD_CK_HPP
#define NEIGHBORHOOD_CK_HPP

#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<>::Neighbor(const ExecutionPolicy &ex_policy,
                     SPHAdaptation *sph_adaptation, DiscreteVariable<Vecd> *dv_pos)
    : kernel_(*sph_adaptation->getKernel()),
      source_pos_(dv_pos->DelegatedDataField(ex_policy)),
      target_pos_(dv_pos->DelegatedDataField(ex_policy)){};
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<>::Neighbor(const ExecutionPolicy &ex_policy,
                     SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
                     DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos)
    : kernel_(*sph_adaptation->getKernel()),
      source_pos_(dv_pos->DelegatedDataField(ex_policy)),
      target_pos_(dv_contact_pos->DelegatedDataField(ex_policy))
{
    KernelWendlandC2CK contact_kernel(*contact_adaptation->getKernel());
    if (kernel_.CutOffRadius() < contact_kernel.CutOffRadius())
    {
        kernel_ = contact_kernel;
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
NeighborList::NeighborList(const ExecutionPolicy &ex_policy,
                           DiscreteVariable<UnsignedInt> *dv_neighbor_index,
                           DiscreteVariable<UnsignedInt> *dv_particle_offset)
    : neighbor_index_(dv_neighbor_index->DelegatedDataField(ex_policy)),
      particle_offset_(dv_particle_offset->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBORHOOD_CK_HPP
