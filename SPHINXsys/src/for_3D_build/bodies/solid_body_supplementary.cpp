#include "solid_body.h"
#include "elastic_solid.h"
#include "solid_particles.h"
#include "sph_system.h"


namespace SPH {
	//===============================================================//
	void SolidBodyPartForSimbody::TagBodyPartParticles()
	{
		SolidBodyPart::TagBodyPartParticles();

		Real body_part_volume(0);
		initial_mass_center_ = Vec3(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData &base_particle_data_i
				= solid_body_->base_particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= dynamic_cast<SolidParticles*>(solid_body_->base_particles_->PointToThisObject())
				->solid_body_data_[index_particle_i];

			initial_mass_center_ += base_particle_data_i.Vol_*solid_data_i.pos_0_;
			body_part_volume += base_particle_data_i.Vol_;
		}

		initial_mass_center_ /= body_part_volume;

		//computing unit intertia
		Vec3 intertia_moments(0);
		Vec3 intertia_products(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData &base_particle_data_i
				= solid_body_->base_particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= dynamic_cast<SolidParticles*>(solid_body_->base_particles_->PointToThisObject())
				->solid_body_data_[index_particle_i];

			Vec3d displacement = (solid_data_i.pos_0_ - initial_mass_center_);
			intertia_moments[0] += base_particle_data_i.Vol_
				*(displacement[1] * displacement[1] + displacement[2] * displacement[2]);
			intertia_moments[1] += base_particle_data_i.Vol_
				*(displacement[0] * displacement[0] + displacement[2] * displacement[2]);
			intertia_moments[2] += base_particle_data_i.Vol_
				*(displacement[0] * displacement[0] + displacement[1] * displacement[1]);
			intertia_products[0] -= base_particle_data_i.Vol_*displacement[0] * displacement[1];
			intertia_products[1] -= base_particle_data_i.Vol_*displacement[0] * displacement[2];
			intertia_products[2] -= base_particle_data_i.Vol_*displacement[1] * displacement[2];

		}
		intertia_moments /= body_part_volume;
		intertia_products /= body_part_volume;

		body_part_mass_properties_
			= new SimTK::MassProperties(body_part_volume*solid_body_density_,
				Vec3(0), UnitInertia(intertia_moments, intertia_products));
	}
	//===============================================================//
}