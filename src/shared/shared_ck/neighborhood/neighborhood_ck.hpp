#ifndef NEIGHBORHOOD_CK_HPP
#define NEIGHBORHOOD_CK_HPP

#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<Base>::Neighbor(const ExecutionPolicy &ex_policy,
                         SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
                         DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos)
    : kernel_(*sph_adaptation->getKernel()),
      kernel_size_square_(pow(kernel_.KernelSize(), 2)),
      source_pos_(dv_pos->DelegatedData(ex_policy)),
      target_pos_(dv_contact_pos->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<>::Neighbor(const ExecutionPolicy &ex_policy,
                     SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
                     DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos)
    : Neighbor<Base>(ex_policy, sph_adaptation, contact_adaptation, dv_pos, dv_contact_pos),
      inv_h_(1.0 / SMAX(sph_adaptation->ReferenceSmoothingLength(),
                        contact_adaptation->ReferenceSmoothingLength())) {}
//=================================================================================================//
template <class SmoothingLengthType>
template <class ExecutionPolicy>
Neighbor<SmoothingLengthType>::Neighbor(
    const ExecutionPolicy &ex_policy,
    SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
    DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos,
    SmoothingLengthType &smoothing_length)
    : Neighbor<Base>(ex_policy, sph_adaptation, contact_adaptation, dv_pos, dv_contact_pos),
      inv_h_(smoothing_length.getInvH()) {}
//=================================================================================================//
template <class SmoothingLengthType>
Neighbor<SmoothingLengthType>::NeighborCriterion::
    NeighborCriterion(Neighbor<SmoothingLengthType> &neighbor)
    : source_pos_(neighbor.source_pos_), target_pos_(neighbor.target_pos_),
      inv_h_(neighbor.inv_h_), kernel_size_square_(neighbor.kernel_size_square_) {}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBORHOOD_CK_HPP
