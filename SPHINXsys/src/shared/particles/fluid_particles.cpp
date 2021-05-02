/** 
 * @file fluid_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "fluid_particles.h"

#include "base_body.h"
#include "weakly_compressible_fluid.h"
#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	FluidParticles::FluidParticles(SPHBody *body, Fluid* fluid)
		: BaseParticles(body, fluid)
	{
		fluid->assignFluidParticles(this);
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
	}
	//=================================================================================================//
	ViscoelasticFluidParticles
		::ViscoelasticFluidParticles(SPHBody *body, Oldroyd_B_Fluid* oldroyd_b_fluid)
		: FluidParticles(body, oldroyd_b_fluid)
	{
		oldroyd_b_fluid->assignViscoelasticFluidParticles(this);
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
}
