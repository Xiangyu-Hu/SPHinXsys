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
		Vecd mass_center = Vecd(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_i = body_part_particles_[i];
			Real particle_volume = solid_particles_->Vol_[index_i];
			mass_center += particle_volume * solid_particles_->pos_0_[index_i];
			body_part_volume += particle_volume;
		}

		mass_center /= body_part_volume;
		initial_mass_center_ = Vec3d(mass_center[0], mass_center[1], 0.0);

		//computing unit intertia
		Real Ix = 0.0;
		Real Iy = 0.0;
		Real Iz = 0.0;
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_i = body_part_particles_[i];
			Vecd particle_position = solid_particles_->pos_0_[index_i];
			Real particle_volume = solid_particles_->Vol_[index_i];

			Real r_x = (particle_position[1] - mass_center[1]);
			Ix += particle_volume * r_x * r_x;
			Real r_y = (particle_position[0] - mass_center[0]);
			Iy += particle_volume * r_y * r_y;
			Iz += particle_volume * (particle_position - mass_center).normSqr();
		}
		Ix /= body_part_volume;
		Iy /= body_part_volume;
		Iz /= body_part_volume;

		body_part_mass_properties_ = mass_properties_ptr_keeper_.createPtr<SimTK::MassProperties>(
			body_part_volume * solid_body_density_, Vec3d(0), SimTK::UnitInertia(Ix, Iy, Iz));
	}
	//=================================================================================================//
}
