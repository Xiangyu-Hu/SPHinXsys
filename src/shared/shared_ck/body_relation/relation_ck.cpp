#include "relation_ck.hpp"

namespace SPH
{
//=================================================================================================//
Relation<Base>::Relation(SPHBody &sph_body)
    : sph_body_(sph_body),
      particles_(sph_body.getBaseParticles()),
      offset_list_size_(particles_.ParticlesBound() + 1) {}
//=================================================================================================//
Relation<Inner<>>::Relation(RealBody &real_body)
    : Relation<Base>(real_body), real_body_(&real_body),
      cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())),
      dv_neighbor_index_(addRelationVariable<UnsignedInt>("NeighborIndex", offset_list_size_)),
      dv_particle_offset_(addRelationVariable<UnsignedInt>("ParticleOffset", offset_list_size_)) {}
//=================================================================================================//
void Relation<Inner<>>::registerComputingKernel(execution::Implementation<Base> *implementation)
{
    all_inner_computing_kernels_.push_back(implementation);
}
//=================================================================================================//
void Relation<Inner<>>::resetComputingKernelUpdated()
{
    for (size_t k = 0; k != all_inner_computing_kernels_.size(); ++k)
    {
        all_inner_computing_kernels_[k]->resetUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
