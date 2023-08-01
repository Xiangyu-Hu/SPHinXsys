#include "porous_solid_particles.h"
#include "porous_media_solid.h"


namespace SPH
{
namespace multi_species_continuum
{	
 
    //=============================================================================================//
	PorousMediaParticles::PorousMediaParticles(SPHBody &body, PorousMediaSolid *porous_solid)
		: ElasticSolidParticles(body, porous_solid), porous_solid_(*porous_solid) {}
	//=============================================================================================//
	void PorousMediaParticles::initializeOtherVariables()
	{
		ElasticSolidParticles::initializeOtherVariables();
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerVariable(fluid_velocity_, "FluidVelocity");//
		registerVariable(relative_fluid_flux_, "RelativeFluidFlux");
 		registerVariable(fluid_saturation_, "FluidSaturation", Real(Eps));
		registerVariable(fluid_mass_, "FluidMass", Real(Eps));

		registerVariable(Vol_update_, "UpdateVolume", [&](size_t i) -> Real { return Vol_[i]; });
		registerVariable(dfluid_mass_dt_, "FluidMassIncrement");
		registerVariable(total_mass_, "TotalMass");
		registerVariable(total_momentum_, "TotalMomentum");
		registerVariable(dtotal_momentum_dt_, "TotalMomentumIncrement");
		
		registerVariable(outer_fluid_velocity_relative_fluid_flux_, "OuterFluidVelocityRelativeFluidFlux");
		registerVariable(Stress_, "Stress");

		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addVariableToWrite<Real>("FluidSaturation");
        addVariableToWrite<Real>("FluidMass");
 
	}
//=================================================================================================//
}
} // namespace SPH
