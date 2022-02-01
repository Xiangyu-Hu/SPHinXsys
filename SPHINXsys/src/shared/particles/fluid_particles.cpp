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
		registerAVariable<indexScalar, Real>(p_, "Pressure");
		registerAVariable<indexScalar, Real>(drho_dt_, "DensityChangeRate");
		registerAVariable<indexScalar, Real>(rho_sum_, "DensitySummation");
		registerAVariable<indexInteger, int>(surface_indicator_, "SurfaceIndicator");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		registerASortableVariable<indexVector, Vecd>("Position");
		registerASortableVariable<indexVector, Vecd>("Velocity");
		registerASortableVariable<indexScalar, Real>("Mass");
		registerASortableVariable<indexScalar, Real>("Density");
		registerASortableVariable<indexScalar, Real>("Pressure");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexScalar, Real>(variables_to_restart_, "Pressure");
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
		registerAVariable<indexMatrix, Matd>(tau_, "ElasticStress");
		registerAVariable<indexMatrix, Matd>(dtau_dt_, "ElasticStressChangeRate");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		registerASortableVariable<indexMatrix, Matd>("ElasticStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexMatrix, Matd>(variables_to_restart_, "ElasticStress");
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
		registerAVariable<indexVector, Vecd>(mom_, "Momentum");
		registerAVariable<indexVector, Vecd>(dmom_dt_, "MomentumChangeRate");
		registerAVariable<indexVector, Vecd>(dmom_dt_prior_, "OtherMomentumChangeRate");
		registerAVariable<indexScalar, Real>(E_, "TotalEnergy");
		registerAVariable<indexScalar, Real>(dE_dt_, "TotalEnergyChangeRate");
		registerAVariable<indexScalar, Real>(dE_dt_prior_, "OtherEnergyChangeRate");
		//----------------------------------------------------------------------
		//		add output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<indexScalar, Real>("Pressure");
	}
	//=================================================================================================//
}
