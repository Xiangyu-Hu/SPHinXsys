/**
 * @file 	solid_body_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_body.h"
#include "solid_particles.h"

namespace SPH {
	//===============================================================//
	void SolidBodyPartForSimbody::TagBodyPartParticles()
	{
		BodyPartByParticle::TagBodyPartParticles();

		Real body_part_volume(0);
		Vecd mass_center = Vecd(0);
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData &base_particle_data_i
				= body_->base_particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= dynamic_cast<SolidParticles*>(body_->base_particles_->PointToThisObject())
				->solid_body_data_[index_particle_i];

			mass_center += base_particle_data_i.Vol_*solid_data_i.pos_0_;
			body_part_volume += base_particle_data_i.Vol_;
		}

		mass_center /= body_part_volume;
		initial_mass_center_ = Vec3(mass_center[0], mass_center[1], 0.0);

		//computing unit intertia
		Real Ix = 0.0;
		Real Iy = 0.0;
		Real Iz = 0.0;
		for (size_t i = 0; i < body_part_particles_.size(); ++i)
		{
			size_t index_particle_i = body_part_particles_[i];
			BaseParticleData &base_particle_data_i
				= body_->base_particles_->base_particle_data_[index_particle_i];
			SolidParticleData &solid_data_i
				= dynamic_cast<SolidParticles*>(body_->base_particles_->PointToThisObject())
				->solid_body_data_[index_particle_i];

			Vecd displacement = (solid_data_i.pos_0_ - mass_center);
			Real r_x = (solid_data_i.pos_0_[1] - mass_center[1]);
			Ix += base_particle_data_i.Vol_*r_x*r_x;
			Real r_y = (solid_data_i.pos_0_[0] - mass_center[0]);
			Iy += base_particle_data_i.Vol_*r_y*r_y;
			Iz += base_particle_data_i.Vol_
				*(solid_data_i.pos_0_ - mass_center).normSqr();
		}
		Ix /= body_part_volume;
		Iy /= body_part_volume;
		Iz /= body_part_volume;

		body_part_mass_properties_ 
			= new SimTK::MassProperties(body_part_volume*solid_body_density_, 
				Vec3(0), UnitInertia(Ix, Iy, Iz));
	}
	//===============================================================//
}