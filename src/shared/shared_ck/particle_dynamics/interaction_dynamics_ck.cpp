#include "interaction_dynamics_ck.h"

namespace SPH
{
//=================================================================================================//
InteractionDynamics<Inner<>>::InteractionDynamics(InnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      sph_adaptation_(sph_body_.sph_adaptation_),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      dv_neighbor_index_(inner_relation.getNeighborIndex()),
      dv_particle_offset_(inner_relation.getParticleOffset()) {}
//=================================================================================================//
InteractionDynamics<Contact<>>::InteractionDynamics(ContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()),
      sph_adaptation_(sph_body_.sph_adaptation_),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      contact_bodies_(contact_relation.getContactBodies()),
      contact_adaptations_(contact_relation.getContactAdaptations()),
      dv_contact_neighbor_index_(contact_relation.getContactNeighborIndex()),
      dv_contact_particle_offset_(contact_relation.getContactParticleOffset()) {}
//=================================================================================================//
} // namespace SPH
