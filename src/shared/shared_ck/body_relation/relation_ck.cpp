#include "relation_ck.hpp"

namespace SPH
{
//=================================================================================================//
Relation<Base>::Relation(SPHBody &sph_body)
    : sph_body_(sph_body),
      particles_(sph_body.getBaseParticles()),
      offset_list_size_(particles_.RealParticlesBound() + 1) {}
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
Relation<Contact<>>::Relation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
    : Relation<Base>(sph_body), contact_bodies_(contact_sph_bodies)
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        const std::string name = contact_bodies_[k]->getName();
        contact_particles_.push_back(&contact_bodies_[k]->getBaseParticles());
        contact_adaptations_.push_back(contact_bodies_[k]->sph_adaptation_);

        CellLinkedList *target_cell_linked_list =
            DynamicCast<CellLinkedList>(this, &contact_bodies_[k]->getCellLinkedList());
        target_cell_linked_lists_.push_back(target_cell_linked_list);

        dv_contact_neighbor_index_.push_back(addRelationVariable<UnsignedInt>(
            "Contact" + name + "NeighborIndex", offset_list_size_));
        dv_contact_particle_offset_.push_back(addRelationVariable<UnsignedInt>(
            "Contact" + name + "ParticleOffset", offset_list_size_));
        all_contact_computing_kernels_.resize(contact_bodies_.size());
    }
}
//=================================================================================================//
void Relation<Contact<>>::registerComputingKernel(
    execution::Implementation<Base> *implementation, UnsignedInt contact_index)
{
    all_contact_computing_kernels_[contact_index].push_back(implementation);
}
//=================================================================================================//
void Relation<Contact<>>::resetComputingKernelUpdated(UnsignedInt contact_index)
{
    for (size_t k = 0; k != all_contact_computing_kernels_[contact_index].size(); ++k)
    {
        all_contact_computing_kernels_[contact_index][k]->resetUpdated();
    }
}
//=================================================================================================//
} // namespace SPH
