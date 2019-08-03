/**
 * @file 	base_particle_generator.h
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The direct generator simply generate
 * 			particle with given position and volume. 
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"

namespace SPH {

	/** Preclaimed classes*/
	class SPHBody;
	class Particles;
	/**
	 * @class ParticleGenerator.
	 * @brief Base abstract class for particle generation.
	 */
	class ParticleGenerator
	{
	public:
		ParticleGenerator(SPHBody &sph_body);
		virtual ~ParticleGenerator() {};
		virtual void CreateBaseParticles(SPHBody &sph_body, Particles &base_particles) = 0;
	};
	/**
	 * @class ParticleGeneratorDirect
	 * @brief Generate particle directly from poistion-and-volume data.
	 */
	class ParticleGeneratorDirect : public ParticleGenerator
	{
	public:
		ParticleGeneratorDirect(SPHBody &sph_body);
		virtual ~ParticleGeneratorDirect() {};
		virtual void CreateBaseParticles(SPHBody &sph_body, Particles &base_particles) override;
	};
}