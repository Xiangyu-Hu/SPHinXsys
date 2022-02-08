/** 
 * @file fluid_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "fluid_particles.h"

#include "base_body.h"
#include "weakly_compressible_fluid.h"
#include "compressible_fluid.h"
#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	FluidParticles::FluidParticles(SPHBody &sph_body,
								   SharedPtr<Fluid> shared_fluid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr)
		: BaseParticles(sph_body, shared_fluid_ptr, particle_generator_ptr)
	{
		shared_fluid_ptr->assignFluidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Real>(p_, "Pressure");
		registerAVariable<Real>(drho_dt_, "DensityChangeRate");
		registerAVariable<Real>(rho_sum_, "DensitySummation");
		registerAVariable<int>(surface_indicator_, "SurfaceIndicator");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		registerASortableVariable<Vecd>("Position");
		registerASortableVariable<Vecd>("Velocity");
		registerASortableVariable<Real>("Mass");
		registerASortableVariable<Real>("Density");
		registerASortableVariable<Real>("Pressure");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<Real>(variables_to_restart_, "Pressure");
	}
	//=================================================================================================//
	ViscoelasticFluidParticles::
		ViscoelasticFluidParticles(SPHBody &sph_body,
								   SharedPtr<Oldroyd_B_Fluid> shared_oldroyd_b_fluid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr)
		: FluidParticles(sph_body, shared_oldroyd_b_fluid_ptr, particle_generator_ptr)
	{
		shared_oldroyd_b_fluid_ptr->assignViscoelasticFluidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Matd>(tau_, "ElasticStress");
		registerAVariable<Matd>(dtau_dt_, "ElasticStressChangeRate");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		registerASortableVariable<Matd>("ElasticStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<Matd>(variables_to_restart_, "ElasticStress");
	}
	//=================================================================================================//
	CompressibleFluidParticles::
		CompressibleFluidParticles(SPHBody &sph_body,
								   SharedPtr<CompressibleFluid> shared_compressiblefluid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr)
		: FluidParticles(sph_body, shared_compressiblefluid_ptr, particle_generator_ptr)
	{
		shared_compressiblefluid_ptr->assignCompressibleFluidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Vecd>(mom_, "Momentum");
		registerAVariable<Vecd>(dmom_dt_, "MomentumChangeRate");
		registerAVariable<Vecd>(dmom_dt_prior_, "OtherMomentumChangeRate");
		registerAVariable<Real>(E_, "TotalEnergy");
		registerAVariable<Real>(dE_dt_, "TotalEnergyChangeRate");
		registerAVariable<Real>(dE_dt_prior_, "OtherEnergyChangeRate");
		//----------------------------------------------------------------------
		//		add output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Real>("Pressure");
	}
	//=================================================================================================//
	WeaklyCompressibleFluidParticles
		::WeaklyCompressibleFluidParticles(SPHBody &sph_body,
								   SharedPtr<Fluid> shared_fluid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr)
		: FluidParticles(sph_body, shared_fluid_ptr, particle_generator_ptr)
	{
		shared_fluid_ptr->assignFluidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Real>(dmass_dt_, "MassChangeRate");
		registerAVariable<Vecd>(mom_, "Momentum");
		registerAVariable<Vecd>(dmom_dt_, "MomentumChangeRate");
		registerAVariable<Vecd>(dmom_dt_prior_, "OtherMomentumChangeRate");
		//----------------------------------------------------------------------
		//		add output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Real>("Pressure");
	}
	//=================================================================================================//
}
