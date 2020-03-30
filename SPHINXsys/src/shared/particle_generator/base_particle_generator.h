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
	class BaseParticles;

	/**
	 * @class ParticleGenerator.
	 * @brief Base abstract class for particle generation.
	 */
	class ParticleGenerator
	{
	protected:
		SPHBody &body_;
	public:
		ParticleGenerator(SPHBody &body);
		virtual ~ParticleGenerator() {};

		/** Create lattice particle for a body. */
		virtual void CreateBaseParticles(BaseParticles* base_particles) = 0;
	};
	/**
	 * @class ParticleGeneratorDirect
	 * @brief Generate particle directly from poistion-and-volume data.
	 */
	class ParticleGeneratorDirect : public ParticleGenerator
	{
	public:
		ParticleGeneratorDirect(SPHBody &body);
		virtual ~ParticleGeneratorDirect() {};
		virtual void CreateBaseParticles(BaseParticles* base_particles) override;
	};
}