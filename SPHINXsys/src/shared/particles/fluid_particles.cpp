#include "fluid_particles.h"

#include "base_body.h"
#include "weakly_compressible_fluid.h"
#include "compressible_fluid.h"
#include "xml_engine.h"
//=====================================================================================================//
namespace SPH
{
	//=================================================================================================//
	FluidParticles::FluidParticles(SPHBody &sph_body, Fluid *fluid)
		: BaseParticles(sph_body, fluid), fluid_(*fluid) {}
	//=================================================================================================//
	void FluidParticles::initializeOtherVariables()
	{
		BaseParticles::initializeOtherVariables();

		registerVariable(p_, "Pressure");
		registerVariable(drho_dt_, "DensityChangeRate");
		registerVariable(rho_sum_, "DensitySummation");
		registerVariable(surface_indicator_, "SurfaceIndicator");
		/**
		 *	register sortable particle data
		 */
		registerSortableVariable<Vecd>("Position");
		registerSortableVariable<Vecd>("Velocity");
		registerSortableVariable<Real>("MassiveMeasure");
		registerSortableVariable<Real>("Density");
		registerSortableVariable<Real>("Pressure");
		registerSortableVariable<Real>("VolumetricMeasure");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addVariableToRestart<Real>("Pressure");
	}
	//=================================================================================================//
	ViscoelasticFluidParticles::
		ViscoelasticFluidParticles(SPHBody &sph_body, Oldroyd_B_Fluid *oldroyd_b_fluid)
		: FluidParticles(sph_body, oldroyd_b_fluid),
		  oldroyd_b_fluid_(*oldroyd_b_fluid) {}
	//=================================================================================================//
	void ViscoelasticFluidParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();

		registerVariable(tau_, "ElasticStress");
		registerVariable(dtau_dt_, "ElasticStressChangeRate");
		/**
		 *	register sortable particle data
		 */		
		registerSortableVariable<Matd>("ElasticStress");
		/**
		 *	add restart output particle data
		 */	
		addVariableToRestart<Matd>("ElasticStress");
	}
	//=================================================================================================//
	CompressibleFluidParticles::
		CompressibleFluidParticles(SPHBody &sph_body, CompressibleFluid *compressible_fluid)
		: FluidParticles(sph_body, compressible_fluid),
		  compressible_fluid_(*compressible_fluid) {}
	//=================================================================================================//
	void CompressibleFluidParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();
		/**
		 *	register sortable particle data
		 */	
		registerVariable(mom_, "Momentum");
		registerVariable(dmom_dt_, "MomentumChangeRate");
		registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
		registerVariable(E_, "TotalEnergy");
		registerVariable(dE_dt_, "TotalEnergyChangeRate");
		registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
		/**
		 *	add output particle data
		 */	
		addVariableToWrite<Real>("Pressure");
	}
	//=================================================================================================//
	WeaklyCompressibleFluidParticles::WeaklyCompressibleFluidParticles(SPHBody &sph_body, Fluid *fluid)
		: FluidParticles(sph_body, fluid) {}
	//=================================================================================================//
	void WeaklyCompressibleFluidParticles::initializeOtherVariables()
	{
		FluidParticles::initializeOtherVariables();
		/**
		 *	register sortable particle data
		 */	
		registerVariable(dmass_dt_, "MassChangeRate");
		registerVariable(mom_, "Momentum");
		registerVariable(dmom_dt_, "MomentumChangeRate");
		registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
		/**
		 *	add output particle data
		 */	
		addVariableToWrite<Real>("Pressure");
	}
	//=================================================================================================//
}
//=====================================================================================================//