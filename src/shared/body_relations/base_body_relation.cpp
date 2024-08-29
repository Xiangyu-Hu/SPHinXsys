#include "base_body_relation.h"

#include "base_particle_dynamics.h"
#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
RealBodyVector BodyPartsToRealBodies(BodyPartVector body_parts)
{
    RealBodyVector real_bodies;

    for (size_t k = 0; k != body_parts.size(); ++k)
    {
        real_bodies.push_back(DynamicCast<RealBody>(body_parts[k], &body_parts[k]->getSPHBody()));
    }
    return real_bodies;
}
//=================================================================================================//
SPHRelation::SPHRelation(SPHBody &sph_body)
    : sph_body_(sph_body),
      base_particles_(sph_body.getBaseParticles()),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position")),
      particle_offset_list_size_(base_particles_.RealParticlesBound() + 1) {}
//=================================================================================================//
BaseInnerRelation::BaseInnerRelation(RealBody &real_body)
    : SPHRelation(real_body), real_body_(&real_body),
      dv_neighbor_index_(addRelationVariable<UnsignedInt>("NeighborIndex", particle_offset_list_size_)),
      dv_particle_offset_(addRelationVariable<UnsignedInt>("ParticleOffset", particle_offset_list_size_))
{
    subscribeToBody();
    inner_configuration_.resize(base_particles_.RealParticlesBound(), Neighborhood());
}
//=================================================================================================//
void BaseInnerRelation::resetNeighborhoodCurrentSize()
{
    parallel_for(
        IndexRange(0, base_particles_.TotalRealParticles()),
        [&](const IndexRange &r)
        {
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                inner_configuration_[num].current_size_ = 0;
            }
        },
        ap);
}
//=================================================================================================//
BaseContactRelation::BaseContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
    : SPHRelation(sph_body), contact_bodies_(contact_sph_bodies)
{
    subscribeToBody();
    contact_configuration_.resize(contact_bodies_.size());
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        const std::string name = contact_bodies_[k]->getName();
        contact_adaptations_.push_back(contact_bodies_[k]->sph_adaptation_);
        contact_configuration_[k].resize(base_particles_.RealParticlesBound(), Neighborhood());
        dv_contact_neighbor_index_.push_back(addRelationVariable<UnsignedInt>(
            "Contact" + name + "NeighborIndex", particle_offset_list_size_));
        dv_contact_particle_offset_.push_back(addRelationVariable<UnsignedInt>(
            "Contact" + name + "ParticleOffset", particle_offset_list_size_));
    }
}
//=================================================================================================//
void BaseContactRelation::resetNeighborhoodCurrentSize()
{
    for (size_t k = 0; k != contact_bodies_.size(); ++k)
    {
        parallel_for(
            IndexRange(0, base_particles_.TotalRealParticles()),
            [&](const IndexRange &r)
            {
                for (size_t num = r.begin(); num != r.end(); ++num)
                {
                    contact_configuration_[k][num].current_size_ = 0;
                }
            },
            ap);
    }
}
//=================================================================================================//
} // namespace SPH
