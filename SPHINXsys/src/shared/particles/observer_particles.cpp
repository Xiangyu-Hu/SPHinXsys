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
	ObserverParticles::ObserverParticles(string body_name)
		: Particles(body_name)
	{

	}
	//===============================================================//
	void ObserverParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		observer_body_data_.push_back(ObserverParticleData());
	}
	//===============================================================//
	ObserverParticles* ObserverParticles::PointToThisObject() 
	{
		return this;
	}
}