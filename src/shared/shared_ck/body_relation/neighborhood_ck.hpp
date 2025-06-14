#ifndef NEIGHBORHOOD_CK_HPP
#define NEIGHBORHOOD_CK_HPP

#include "neighborhood_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
Neighbor<Base>::Neighbor(
    const ExecutionPolicy &ex_policy,
    SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
    DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos)
    : kernel_(*sph_adaptation->getKernel()),
      source_pos_(dv_pos->DelegatedData(ex_policy)),
      target_pos_(dv_contact_pos->DelegatedData(ex_policy))
{
    KernelTabulatedCK contact_kernel(*contact_adaptation->getKernel());
    if (kernel_.CutOffRadius() < contact_kernel.CutOffRadius())
    {
        kernel_ = contact_kernel;
    }
}
//=================================================================================================//
template <class NeighborMethod>
template <class ExecutionPolicy>
Neighbor<NeighborMethod>::Neighbor(
    const ExecutionPolicy &ex_policy,
    SPHAdaptation *sph_adaptation, SPHAdaptation *contact_adaptation,
    DiscreteVariable<Vecd> *dv_pos, DiscreteVariable<Vecd> *dv_contact_pos,
    NeighborMethod &neighbor_method)
    : Neighbor<Base>(ex_policy, sph_adaptation, contact_adaptation, dv_pos, dv_contact_pos),
      cut_radius_square_(kernel_.CutOffRadiusSqr()) {}
//=================================================================================================//
} // namespace SPH
#endif // NEIGHBORHOOD_CK_HPP
