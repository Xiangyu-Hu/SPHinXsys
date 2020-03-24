/**
 * @file observer_particles.cpp
 * @author	Chi Zhang and Xiangyu Hu
 * @version	0.1
 */
#include "observer_particles.h"
#include "base_material.h"

namespace SPH
{
	//===============================================================//
	ObserverParticleData::ObserverParticleData()
	{

	}
	//===============================================================//
	ObserverParticles::ObserverParticles(SPHBody *body)
		: BaseParticles(body, new BaseMaterial())
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