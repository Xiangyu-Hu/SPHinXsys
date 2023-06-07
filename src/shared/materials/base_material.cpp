#include "base_material.h"
#include "base_particles.hpp"

namespace SPH
{
	//=================================================================================================//
	void BaseMaterial::registerReloadLocalParameters(BaseParticles* base_particles)
	{
		base_particles->registerVariable(surface_indicator_, "SurfaceIndicator");
	}
	
	void Fluid::registerReloadLocalParameters(BaseParticles* base_particles)
	{
		BaseMaterial::registerReloadLocalParameters(base_particles);

		base_particles->registerVariable(p_, "Pressure");
		base_particles->registerVariable(drho_dt_, "DensityChangeRate");
		base_particles->registerVariable(rho_sum_, "DensitySummation");
	}
	//=================================================================================================//
}
