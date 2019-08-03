/**
 * @file weakly_compressible_fluid_particles.cpp
 * @brief Definition of funcitons declared in weakly_compressible_fluid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "fluid_particles.h"

namespace SPH
{
	//===============================================================//
	FluidParticleData::FluidParticleData()
		: mass_(0.0), sigma_0_(0.0), rho_0_(1.0), rho_n_(1.0), p_(0.0), 
		drho_dt_(0.0), div_correction_(1.0), vel_trans_(0),
		dvel_dt_trans_(0), c_(1.0), vorticity_(0), vort_2d_(0.0)
	{

	}
	//===============================================================//
	FluidParticles::FluidParticles(SPHBody *body) : Particles(body)
	{
		for(size_t i = 0; i < base_particle_data_.size(); ++i)
			fluid_particle_data_.push_back(FluidParticleData());
	}
	//===============================================================//
	FluidParticles* FluidParticles::PointToThisObject() 
	{
		return this;
	}
	//===============================================================//
	ViscoelasticFluidParticleData::ViscoelasticFluidParticleData()
		: tau_(0), dtau_dt_(0)
	{

	}
	//===============================================================//
	ViscoelasticFluidParticles::ViscoelasticFluidParticles(SPHBody *body)
		: FluidParticles(body)
	{
		for (size_t i = 0; i < base_particle_data_.size(); ++i)
			viscoelastic_particle_data_.push_back(ViscoelasticFluidParticleData());
	}
	//===============================================================//
	ViscoelasticFluidParticles* ViscoelasticFluidParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
}