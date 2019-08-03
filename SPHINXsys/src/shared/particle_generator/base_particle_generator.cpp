/**
 * @file 	base_particle_generator.cpp
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The direct generator simply generate
 * 			particle with given position and volume. The lattice generator generate
 * 			at lattice position by check whether the poision is contained by a SPH body.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#include "base_particle_generator.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//---------------------------------------------------------------//
	ParticleGenerator
		::ParticleGenerator(SPHBody &sph_body)
	{

	}
	//---------------------------------------------------------------//
	ParticleGeneratorDirect
		::ParticleGeneratorDirect(SPHBody &sph_body)
		: ParticleGenerator(sph_body)
	{

	}
	//---------------------------------------------------------------//
	void ParticleGeneratorDirect
		::CreateBaseParticles(SPHBody &sph_body, Particles &base_particles)
	{
		size_t number_of_particles = 0;
		for (int i = 0; i < sph_body.body_input_points_volumes_.size(); ++i) 
		{
			base_particles.InitializeABaseParticle(sph_body.body_input_points_volumes_[i].first,
				sph_body.body_input_points_volumes_[i].second);
			number_of_particles++;
		}

		sph_body.number_of_particles_ = number_of_particles;
	}
	//---------------------------------------------------------------//
}