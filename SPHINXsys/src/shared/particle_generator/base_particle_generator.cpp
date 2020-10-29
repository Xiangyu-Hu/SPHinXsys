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
	void ParticleGenerator::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
	}
	//=================================================================================================//
	void ParticleGeneratorDirect
		::CreateBaseParticles(BaseParticles* base_particles)
	{
		size_t number_of_particles = 0;
		auto& body_input_points_volumes = sph_body_->body_input_points_volumes_;
		for (size_t i = 0; i < body_input_points_volumes.size(); ++i)
		{
			base_particles->initializeABaseParticle(body_input_points_volumes[i].first,
				body_input_points_volumes[i].second);
			number_of_particles++;
		}

		sph_body_->number_of_particles_ = number_of_particles;
	}
	//=================================================================================================//
}
