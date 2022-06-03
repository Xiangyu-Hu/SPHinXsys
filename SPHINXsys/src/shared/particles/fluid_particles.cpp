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
	FluidParticles::FluidParticles(SPHBody &sph_body, Fluid *fluid)
		: BaseParticles(sph_body, fluid) {}
	//=================================================================================================//
	void FluidParticles::initializeOtherVariables()
	{
		BaseParticles::initializeOtherVariables();

		registerAVariable(p_, "Pressure");
		registerAVariable(drho_dt_, "DensityChangeRate");
		registerAVariable(rho_sum_, "DensitySummation");
		registerAVariable(surface_indicator_, "SurfaceIndicator");
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
		addAVariableToRestart<Real>("Pressure");
	}
	//=================================================================================================//
	ViscoelasticFluidParticles::
		ViscoelasticFluidParticles(SPHBody &sph_body, Oldroyd_B_Fluid *oldroyd_b_fluid)
		: FluidParticles(sph_body, oldroyd_b_fluid) {}
	//=================================================================================================//
	void ViscoelasticFluidParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();

		registerAVariable(tau_, "ElasticStress");
		registerAVariable(dtau_dt_, "ElasticStressChangeRate");
		//----------------------------------------------------------------------
		//		register sortable particle data
		//----------------------------------------------------------------------
		registerASortableVariable<Matd>("ElasticStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableToRestart<Matd>("ElasticStress");
	}
	//=================================================================================================//
	CompressibleFluidParticles::
		CompressibleFluidParticles(SPHBody &sph_body, CompressibleFluid *compressible_fluid)
		: FluidParticles(sph_body, compressible_fluid) {}
	//=================================================================================================//
	void CompressibleFluidParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();

		registerAVariable(mom_, "Momentum");
		registerAVariable(dmom_dt_, "MomentumChangeRate");
		registerAVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
		registerAVariable(E_, "TotalEnergy");
		registerAVariable(dE_dt_, "TotalEnergyChangeRate");
		registerAVariable(dE_dt_prior_, "OtherEnergyChangeRate");
		//----------------------------------------------------------------------
		//		add output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Real>("Pressure");
	}
	//=================================================================================================//
	WeaklyCompressibleFluidParticles::
		WeaklyCompressibleFluidParticles(SPHBody &sph_body, Fluid *fluid)
		: FluidParticles(sph_body, fluid) {}
	//=================================================================================================//
	void WeaklyCompressibleFluidParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();

		registerAVariable(dmass_dt_, "MassChangeRate");
		registerAVariable(mom_, "Momentum");
		registerAVariable(dmom_dt_, "MomentumChangeRate");
		registerAVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
		//----------------------------------------------------------------------
		//		add output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Real>("Pressure");
	}
	//=================================================================================================//
}
