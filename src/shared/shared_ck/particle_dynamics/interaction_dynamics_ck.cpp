#include "interaction_dynamics_ck.h"

namespace SPH
{
//=================================================================================================//
InteractionDynamics<Inner<>>::InteractionDynamics(InnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      kernel_(inner_relation.getSmoothingKernel()),
      dv_neighbor_index(inner_relation.getNeighborIndex()),
      dv_particle_offset(inner_relation.getParticleOffset()) {}
//=================================================================================================//
InteractionDynamics<Contact<>>::InteractionDynamics(ContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")),
      kernel_(contact_relation.getContactSmoothingKernel()),
      dv_neighbor_index(contact_relation.getContactNeighborIndex()),
      dv_particle_offset(contact_relation.getContactParticleOffset()) {}
//=================================================================================================//
} // namespace SPH
