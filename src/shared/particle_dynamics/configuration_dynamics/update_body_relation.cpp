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
      dv_particle_offset_(inner_relation.getParticleOffset()),
      dv_neighbor_size_(DiscreteVariable<UnsignedInt>("NeighborSize", particles_->RealParticlesBound()))
{
    particles_->addVariableToWrite(&dv_neighbor_size_);
}
//=================================================================================================//
} // namespace SPH
