#include "relation_ck.hpp"

namespace SPH
{
//=================================================================================================//
Relation<Inner<>>::Relation(RealBody &real_body)
    : Relation<Base>(real_body, StdVec<RealBody *>{&real_body}),
      real_body_(&real_body), dv_neighbor_index_(dv_target_neighbor_index_[0]),
      dv_particle_offset_(dv_target_particle_offset_[0]) {}
//=================================================================================================//
void Relation<Inner<>>::registerComputingKernel(execution::Implementation<Base> *implementation)
{
    registered_computing_kernels_[0].push_back(implementation);
}
//=================================================================================================//
void Relation<Inner<>>::resetComputingKernelUpdated()
{
    auto &all_inner_computing_kernels = registered_computing_kernels_[0];
    for (size_t k = 0; k != all_inner_computing_kernels.size(); ++k)
    {
        all_inner_computing_kernels[k]->resetUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
