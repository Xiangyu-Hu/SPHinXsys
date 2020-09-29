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
 * @file    fluid_body.h
 * @brief 	This is the class for bodies used for fluid.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_body.h"

#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @brief preclaimed class.
	 */
	class SPHSystem;
	/**
	 * @class FluidBody
	 * @brief Declaration of fluid body.
	 */
	class FluidBody : public RealBody
	{
	public:
		explicit FluidBody(SPHSystem &system, string body_name, int refinement_level,
			ParticleGenerator* particle_generator = new ParticleGeneratorLattice());
		virtual ~FluidBody() {};

		/** The pointer to derived class object. */
		virtual FluidBody* pointToThisObject() override { return this; };
	};
}
