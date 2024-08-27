#include "update_body_relation.h"

namespace SPH
{
//=================================================================================================//
BodyRelationUpdate<Inner<>>::BodyRelationUpdate(InnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      cell_linked_list_(inner_relation.getCellLinkedList()),
      particle_offset_list_size_(inner_relation.getParticleOffsetListSize()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_neighbor_index_(inner_relation.getNeighborIndex()),
      dv_particle_offset_(inner_relation.getParticleOffset())
{
    particles_->addVariableToWrite(dv_particle_offset_);
}
//=================================================================================================//
BodyRelationUpdate<Contact<>>::BodyRelationUpdate(ContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()),
      particle_offset_list_size_(contact_relation.getParticleOffsetListSize()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      contact_cell_linked_list_(contact_relation.getContactCellLinkedList()),
      dv_contact_neighbor_index_(contact_relation.getContactNeighborIndex()),
      dv_contact_particle_offset_(contact_relation.getContactParticleOffset())
{
    for (size_t k = 0; k != dv_contact_particle_offset_.size(); ++k)
    {
        particles_->addVariableToWrite(dv_contact_particle_offset_[k]);
    }
}
//=================================================================================================//
} // namespace SPH
