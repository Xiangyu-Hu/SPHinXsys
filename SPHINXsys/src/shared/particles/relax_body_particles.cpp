/**
 * @file relax_body_particles.cpp
 * @brief Definition of funcitons declared in relax_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "relax_body_particles.h"

namespace SPH {
	//===============================================================//
	RelaxBodyParticleData::RelaxBodyParticleData(Vecd position)
	:pos_0_(position)
	{

	}
	//===============================================================//
	RelaxBodyParticles::RelaxBodyParticles(string body_name)
		: Particles(body_name)
	{

	}
	//===============================================================//
	void RelaxBodyParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		relax_body_data_.push_back(RelaxBodyParticleData(pnt));
	}
	//===============================================================//
	RelaxBodyParticles* RelaxBodyParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
}
