#include "update_body_relation.h"

namespace SPH
{
//=================================================================================================//
BodyRelationUpdate<Inner<>>::BodyRelationUpdate(InnerRelation &inner_relation)
    : LocalInteractionDynamics<Inner<>>(inner_relation),
      cell_linked_list_(inner_relation.getCellLinkedList()),
      particle_offset_list_size_(inner_relation.getParticleOffsetListSize())
{
    particles_->addVariableToWrite(dv_particle_offset_);
}
//=================================================================================================//
BodyRelationUpdate<Contact<>>::BodyRelationUpdate(ContactRelation &contact_relation)
    : LocalInteractionDynamics<Contact<>>(contact_relation),
      particle_offset_list_size_(contact_relation.getParticleOffsetListSize()),
      contact_cell_linked_list_(contact_relation.getContactCellLinkedList())
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        particles_->addVariableToWrite(dv_contact_particle_offset_[k]);
    }
}
//=================================================================================================//
} // namespace SPH
