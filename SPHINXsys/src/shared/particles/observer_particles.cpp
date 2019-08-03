/**
 * @file observer_particles.cpp
 * @brief Definition of funcitons declared in observer_particles.h
 * @author	Chi Zhang and Xiangyu Hu
 * @version	0.1
 */
#include "observer_particles.h"

namespace SPH
{
	//===============================================================//
	ObserverParticleData::ObserverParticleData()
	{

	}
	//===============================================================//
	ObserverParticles::ObserverParticles(SPHBody *body)
		: Particles(body)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
			observer_body_data_.push_back(ObserverParticleData());
	}
	//===============================================================//
	ObserverParticles* ObserverParticles::PointToThisObject() 
	{
		return this;
	}
}