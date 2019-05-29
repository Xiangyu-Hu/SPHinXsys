/**
 * @file solid_body_particles.cpp
 * @brief Definition of funcitons declared in solid_bdoy_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "solid_body_particles.h"

namespace SPH {
	//===============================================================//
	SolidBodyParticleData::SolidBodyParticleData(Vecd position)
		: pos_0_(position), n_0_(0), n_(0), vel_ave_(0), dvel_dt_ave_(0),
		viscous_force_from_fluid_(0), force_from_fluid_(0)
	{

	}
	//===============================================================//
	ElasticBodyParticleData::ElasticBodyParticleData()
		: rho_0_(1.0), rho_n_(1.0), F_(1.0), dF_dt_(0), stress_(0), mass_(1.0),
		B_(0), local_G_(1.0), local_lambda_(1.0), local_eta_(0.0), 
		local_c_(1.0), pos_temp_(0)
	{

	}
	//===============================================================//
	MuscleBodyData
		::MuscleBodyData()
	{

	}
	//===============================================================//
	SolidBodyParticles::SolidBodyParticles(string body_name)
		: Particles(body_name)
	{

	}
	//===============================================================//
	void SolidBodyParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		solid_body_data_.push_back(SolidBodyParticleData(pnt));
	}
	//===============================================================//
	SolidBodyParticles* SolidBodyParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
	ElasticBodyParticles::ElasticBodyParticles(string body_name)
		: SolidBodyParticles(body_name)
	{

	}
	//===============================================================//
	void ElasticBodyParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		solid_body_data_.push_back(SolidBodyParticleData(pnt));
		elastic_body_data_.push_back(ElasticBodyParticleData());
	}
	//===============================================================//
	ElasticBodyParticles* ElasticBodyParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
	MuscleBodyParticles
		::MuscleBodyParticles(string body_name)
		: ElasticBodyParticles(body_name)
	{

	}
	//===============================================================//
	void MuscleBodyParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		solid_body_data_.push_back(SolidBodyParticleData(pnt));
		elastic_body_data_.push_back(ElasticBodyParticleData());
		muscle_body_data_.push_back(MuscleBodyData());
	}
	//===============================================================//
	MuscleBodyParticles* MuscleBodyParticles::PointToThisObject()
	{
		return this;
	}
}
