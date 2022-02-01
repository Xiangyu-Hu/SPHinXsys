/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	compressible_fluid.h
 * @brief 	Describe the compressible fluid which is used 
 * 			model compressible fluids. Here, we have ideal gas equation of states. 
 * @author  Xiangyu Hu, Luhui Han, Chi Zhang and Zhentong Wang
 */

#pragma once

#include "base_material.h"

namespace SPH
{

	class CompressibleFluidParticles;

	/**
	 * @class CompressibleFluid
	 * @brief Ideal gas equation of state (EOS).
	 */
	class CompressibleFluid : public Fluid
	{
	protected:
		Real gamma_; /** heat capacity ratio */
		CompressibleFluidParticles *compressible_fluid_particles_;

	public:
		explicit CompressibleFluid(Real rho0, Real gamma, Real mu = 0.0)
			: Fluid(rho0, mu), gamma_(gamma), compressible_fluid_particles_(nullptr)
		{
			material_type_ = "CompressibleFluid";
		};
		virtual ~CompressibleFluid(){};

		void assignCompressibleFluidParticles(CompressibleFluidParticles *compressible_fluid_particles)
		{
			compressible_fluid_particles_ = compressible_fluid_particles;
		};
		Real HeatCapacityRatio() { return gamma_; };
		virtual Real getPressure(Real rho, Real rho_e) override;
		virtual Real getPressure(Real rho) override { return 0.0; };
		virtual Real DensityFromPressure(Real p) override { return 0.0; };
		virtual Real getSoundSpeed(Real p, Real rho) override;
		virtual CompressibleFluid *ThisObjectPtr() override { return this; };
	};
}
