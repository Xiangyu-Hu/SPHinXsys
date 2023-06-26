#include "base_body_relation.h"
#include "base_particle_dynamics.h"

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
		: sph_body_(sph_body), base_particles_(sph_body.getBaseParticles()) {}
	//=================================================================================================//
	BaseInnerRelation::BaseInnerRelation(RealBody &real_body)
		: SPHRelation(real_body), real_body_(&real_body)
	{
		subscribeToBody();
		resizeConfiguration();
	}
	//=================================================================================================//
	void BaseInnerRelation::resizeConfiguration()
	{
		size_t updated_size = base_particles_.real_particles_bound_;
		inner_configuration_.resize(updated_size, Neighborhood());
	}
	//=================================================================================================//
	void BaseInnerRelation::resetNeighborhoodCurrentSize()
	{
		parallel_for(
			IndexRange(0, base_particles_.total_real_particles_),
			[&](const IndexRange &r)
			{
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					inner_configuration_[num].current_size_ = 0;
				}
			},
			ap);
	}

    void BaseInnerRelation::allocateInnerConfigurationDevice() {
        inner_configuration_device_ = makeSharedDevice<StdSharedVec<NeighborhoodDevice>>(inner_configuration_.size(),
                execution::executionQueue.getQueue());
    }

    void BaseInnerRelation::copyInnerConfigurationToDevice() {
        for (std::size_t i = 0; i < inner_configuration_.size(); ++i)
            inner_configuration_device_->at(i).copyFrom(inner_configuration_.at(i));
    }

    void BaseInnerRelation::copyInnerConfigurationFromDevice() {
        for (std::size_t i = 0; i < inner_configuration_.size(); ++i)
            inner_configuration_device_->at(i).copyTo(inner_configuration_.at(i));
    }

    //=================================================================================================//
	BaseContactRelation::BaseContactRelation(SPHBody &sph_body, RealBodyVector contact_sph_bodies)
		: SPHRelation(sph_body), contact_bodies_(contact_sph_bodies)
	{
		subscribeToBody();
		contact_configuration_.resize(contact_bodies_.size());
	}
	//=================================================================================================//
	void BaseContactRelation::resizeConfiguration()
	{
		size_t updated_size = base_particles_.real_particles_bound_;
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			contact_configuration_[k].resize(updated_size, Neighborhood());
		}
	}
	//=================================================================================================//
	void BaseContactRelation::resetNeighborhoodCurrentSize()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			parallel_for(
				IndexRange(0, base_particles_.total_real_particles_),
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
}
