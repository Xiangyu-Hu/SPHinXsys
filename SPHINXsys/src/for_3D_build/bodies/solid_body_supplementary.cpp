/**
 * @file 	base_body_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
	//=================================================================================================//
	void SolidBodyPartForSimbody::setMassProperties()
	{
		Real body_part_volume(0);
		initial_mass_center_ = Vec3d(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_i = body_part_particles_[i];
			Vecd particle_position = solid_particles_->pos_0_[index_i];
			Real particle_volume = solid_particles_->Vol_[index_i];

			initial_mass_center_ += particle_volume * particle_position;
			body_part_volume += particle_volume;
		}

		initial_mass_center_ /= body_part_volume;

		//computing unit intertia
		Vec3d intertia_moments(0);
		Vec3d intertia_products(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_i = body_part_particles_[i];
			Vecd particle_position = solid_particles_->pos_0_[index_i];
			Real particle_volume = solid_particles_->Vol_[index_i];

			Vec3d displacement = (particle_position - initial_mass_center_);
			intertia_moments[0] += particle_volume * (displacement[1] * displacement[1] + displacement[2] * displacement[2]);
			intertia_moments[1] += particle_volume * (displacement[0] * displacement[0] + displacement[2] * displacement[2]);
			intertia_moments[2] += particle_volume * (displacement[0] * displacement[0] + displacement[1] * displacement[1]);
			intertia_products[0] -= particle_volume * displacement[0] * displacement[1];
			intertia_products[1] -= particle_volume * displacement[0] * displacement[2];
			intertia_products[2] -= particle_volume * displacement[1] * displacement[2];
		}
		intertia_moments /= body_part_volume;
		intertia_products /= body_part_volume;

		body_part_mass_properties_ = mass_properties_ptr_keeper_.createPtr<SimTK::MassProperties>(
			body_part_volume * solid_body_density_, Vec3d(0), SimTK::UnitInertia(intertia_moments, intertia_products));
	}
	//=================================================================================================//
}
