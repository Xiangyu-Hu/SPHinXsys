#include "base_body_relation.h"
#include "base_particle_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	SPHBodyRelation::SPHBodyRelation(SPHBody &sph_body)
		: sph_body_(sph_body), base_particles_(sph_body.base_particles_) {}
	//=================================================================================================//
	BaseBodyRelationInner::BaseBodyRelationInner(RealBody &real_body)
		: SPHBodyRelation(real_body), real_body_(&real_body)
	{
		subscribeToBody();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void BaseBodyRelationInner::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_.base_particles_->real_particles_bound_;
		inner_configuration_.resize(updated_size, Neighborhood());
	}
	//=================================================================================================//
	void BaseBodyRelationInner::resetNeighborhoodCurrentSize()
	{
		parallel_for(
			blocked_range<size_t>(0, base_particles_->total_real_particles_),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					inner_configuration_[num].current_size_ = 0;
				}
			},
			ap);
	}
	//=================================================================================================//
	BaseBodyRelationContact::BaseBodyRelationContact(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
		: SPHBodyRelation(sph_body), contact_bodies_(contact_sph_bodies)
	{
		subscribeToBody();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	BaseBodyRelationContact::BaseBodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts)
		: SPHBodyRelation(sph_body)
	{
		for (size_t k = 0; k != contact_body_parts.size(); ++k)
		{
			contact_bodies_.push_back(DynamicCast<RealBody>(this, &contact_body_parts[k]->getSPHBody()));
		}
		subscribeToBody();
		updateConfigurationMemories();
	}
	//=================================================================================================//
	void BaseBodyRelationContact::updateConfigurationMemories()
	{
		size_t updated_size = sph_body_.base_particles_->real_particles_bound_;
		contact_configuration_.resize(contact_bodies_.size());
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			contact_configuration_[k].resize(updated_size, Neighborhood());
		}
	}
	//=================================================================================================//
	void BaseBodyRelationContact::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			parallel_for(
				blocked_range<size_t>(0, base_particles_->total_real_particles_),
				[&](const blocked_range<size_t> &r)
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
}
