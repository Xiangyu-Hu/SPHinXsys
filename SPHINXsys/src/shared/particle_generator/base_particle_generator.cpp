/**
 * @file 	base_particle_generator.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "base_particle_generator.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//=================================================================================================//
	ParticleGenerator
		::ParticleGenerator(SPHBody &body)
		: body_(body)
	{

	}
	//=================================================================================================//
	ParticleGeneratorDirect
		::ParticleGeneratorDirect(SPHBody &body)
		: ParticleGenerator(body)
	{

	}
	//=================================================================================================//
	void ParticleGeneratorDirect
		::CreateBaseParticles(BaseParticles* base_particles)
	{
		size_t number_of_particles = 0;
		for (size_t i = 0; i < body_.body_input_points_volumes_.size(); ++i) 
		{
			base_particles->InitializeABaseParticle(body_.body_input_points_volumes_[i].first,
				body_.body_input_points_volumes_[i].second, 1.0 / body_.body_input_points_volumes_[i].second);
			number_of_particles++;
		}

		body_.number_of_particles_ = number_of_particles;
	}
	//=================================================================================================//
}
