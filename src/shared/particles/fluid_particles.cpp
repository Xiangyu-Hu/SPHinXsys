 

#include "base_body.h"
#include "weakly_compressible_fluid.h"
#include "compressible_fluid.h"
#include "xml_engine.h"
//=====================================================================================================//
namespace SPH
{
	//=================================================================================================//
	// BaseParticles::BaseParticles(SPHBody &sph_body, Fluid *fluid)
	// 	: BaseParticles(sph_body, fluid), fluid_(*fluid) {}
	// //=================================================================================================//
	// void BaseParticles::initializeOtherVariables()
	// {
	// 	BaseParticles::initializeOtherVariables();

	// 	registerVariable(p_, "Pressure");
	// 	registerVariable(drho_dt_, "DensityChangeRate");
	// 	registerVariable(rho_sum_, "DensitySummation");
	// 	registerVariable(surface_indicator_, "SurfaceIndicator");
	// 	/**
	// 	 *	register sortable particle data
	// 	 */
	// 	registerSortableVariable<Vecd>("Position");
	// 	registerSortableVariable<Vecd>("Velocity");
	// 	registerSortableVariable<Real>("MassiveMeasure");
	// 	registerSortableVariable<Real>("Density");
	// 	registerSortableVariable<Real>("Pressure");
	// 	registerSortableVariable<Real>("VolumetricMeasure");
	// 	//----------------------------------------------------------------------
	// 	//		add restart output particle data
	// 	//----------------------------------------------------------------------
	// 	addVariableToRestart<Real>("Pressure");
	// }
	// //=================================================================================================//
	// ViscoelasticBaseParticles::
	// 	ViscoelasticBaseParticles(SPHBody &sph_body, Oldroyd_B_Fluid *oldroyd_b_fluid)
	// 	: BaseParticles(sph_body, oldroyd_b_fluid),
	// 	  oldroyd_b_fluid_(*oldroyd_b_fluid) {}
	// //=================================================================================================//
	// void ViscoelasticBaseParticles::initializeOtherVariables()
	// {
	// 	BaseParticles::initializeOtherVariables();

	// 	registerVariable(tau_, "ElasticStress");
	// 	registerVariable(dtau_dt_, "ElasticStressChangeRate");
	// 	/**
	// 	 *	register sortable particle data
	// 	 */		
	// 	registerSortableVariable<Matd>("ElasticStress");
	// 	/**
	// 	 *	add restart output particle data
	// 	 */	
	// 	addVariableToRestart<Matd>("ElasticStress");
	// }
	// //=================================================================================================//
	// CompressibleBaseParticles::
	// 	CompressibleBaseParticles(SPHBody &sph_body, CompressibleFluid *compressible_fluid)
	// 	: BaseParticles(sph_body, compressible_fluid),
	// 	  compressible_fluid_(*compressible_fluid) {}
	// //=================================================================================================//
	// void CompressibleBaseParticles::initializeOtherVariables()
	// {
	// 	BaseParticles::initializeOtherVariables();
	// 	/**
	// 	 *	register sortable particle data
	// 	 */	
	// 	registerVariable(mom_, "Momentum");
	// 	registerVariable(dmom_dt_, "MomentumChangeRate");
	// 	registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
	// 	registerVariable(E_, "TotalEnergy");
	// 	registerVariable(dE_dt_, "TotalEnergyChangeRate");
	// 	registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
	// 	/**
	// 	 *	add output particle data
	// 	 */	
	// 	addVariableToWrite<Real>("Pressure");
	// }
	//=================================================================================================//
}
//=====================================================================================================//