#ifndef NEIGHBORHOOD_CK_HPP
#define NEIGHBORHOOD_CK_HPP

#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<>::Neighbor(const ExecutionPolicy &ex_policy,
                     SPHAdaptation *sph_adaptation, DiscreteVariable<Vecd> *dv_pos)
    : source_pos_(dv_pos->DelegatedDataField(ex_policy)),
      target_pos_(dv_pos->DelegatedDataField(ex_policy)),
      kernel_(*sph_adaptation->getKernel()){};
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<>::Neighbor(const ExecutionPolicy &ex_policy,
                     SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
                     DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos)
    : source_pos_(dv_pos->DelegatedDataField(ex_policy)),
      target_pos_(dv_contact_pos->DelegatedDataField(ex_policy))
{
    Kernel *kernel = sph_adaptation->getKernel();
    Kernel *contact_kernel = contact_adaptation->getKernel();
    kernel_ = kernel->SmoothingLength() > contact_kernel->SmoothingLength()
                  ? KernelWendlandC2CK(*kernel)
                  : KernelWendlandC2CK(*contact_kernel);
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
