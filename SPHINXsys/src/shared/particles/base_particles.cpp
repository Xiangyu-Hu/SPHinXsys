/**
 * @file base_particles.cpp
 * @brief Definition of funcitons declared in base_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "base_particles.h"

namespace SPH
{
	//===============================================================//
	int BaseParticleData::ID_max_ = 0;
	//===============================================================//
	BaseParticleData::BaseParticleData(Vecd position, Real volume) 
		: vel_n_(0), pos_n_(position), Vol_(volume), Vol_0_(volume),
		cell_location_(0), dvel_dt_others_(0), dvel_dt_(0)
	{
		particle_ID_ = ID_max_;
		ID_max_++;
	}
	//===============================================================//
	Particles::Particles(string body_name) 
		: body_name_(body_name)
	{

	}
	//===============================================================//
	Particles* Particles::PointToThisObject()
	{
		return this;
	}
}