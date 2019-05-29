/**
 * @file weakly_compressible_fluid_particles.cpp
 * @brief Definition of funcitons declared in weakly_compressible_fluid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "weakly_compressible_fluid_particles.h"

namespace SPH
{
	//===============================================================//
	WeaklyCompressibleFluidParticleData::WeaklyCompressibleFluidParticleData()
		: mass_(0.0), sigma_0_(0.0), rho_0_(1.0), rho_n_(1.0), p_(0.0), 
		drho_dt_(0.0), div_correction_(1.0), vel_trans_(0),
		dvel_dt_trans_(0), c_(1.0), vorticity_(0), vort_2d_(0.0)
	{

	}
	//===============================================================//
	WeaklyCompressibleFluidParticles::WeaklyCompressibleFluidParticles(string body_name)
		: Particles(body_name)
	{

	}
	//===============================================================//
	void WeaklyCompressibleFluidParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		fluid_data_.push_back(WeaklyCompressibleFluidParticleData());
	}
	//===============================================================//
	WeaklyCompressibleFluidParticles* WeaklyCompressibleFluidParticles::PointToThisObject() 
	{
		return this;
	}
	//===============================================================//
	Oldroyd_B_FluidParticleData::Oldroyd_B_FluidParticleData()
		: tau_(0), dtau_dt_(0)
	{

	}
	//===============================================================//
	Oldroyd_B_FluidParticles::Oldroyd_B_FluidParticles(string body_name)
		: WeaklyCompressibleFluidParticles(body_name)
	{

	}
	//===============================================================//
	void Oldroyd_B_FluidParticles::InitializeAParticle(Vecd pnt, Real particle_volume)
	{
		base_particle_data_.push_back(BaseParticleData(pnt, particle_volume));
		fluid_data_.push_back(WeaklyCompressibleFluidParticleData());
		oldroyd_b_data_.push_back(Oldroyd_B_FluidParticleData());
	}
	//===============================================================//
	Oldroyd_B_FluidParticles* Oldroyd_B_FluidParticles::PointToThisObject()
	{
		return this;
	}
	//===============================================================//
}