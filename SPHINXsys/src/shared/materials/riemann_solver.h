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
 * @file 	riemann_solvers.h
 * @brief 	This is the collection of Riemann solvers.
 * @author	Xiangyu Hu
 */

#pragma once

#include "base_data_package.h"

namespace SPH
{
	struct FluidState
	{
		Vecd vel_;
		Real rho_, p_;
		FluidState(Real& rho, Vecd& vel, Real& p) :
			vel_(vel), rho_(rho), p_(p) {};
	};

	class Fluid;

	class NoRiemannSolver
	{
		Fluid& fluid_l_, & fluid_r_;
	public:
		NoRiemannSolver(Fluid& fluid_i, Fluid& fluid_j) :
			fluid_l_(fluid_i), fluid_r_(fluid_j) {};
		Real getPStar(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
		Vecd getVStar(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
	};

	class AcousticRiemannSolver
	{
		Fluid& fluid_i_, & fluid_j_;
	public:
		AcousticRiemannSolver(Fluid& fluid_i, Fluid& fluid_j) :
			fluid_i_(fluid_i), fluid_j_(fluid_j) {};
		Real getPStar(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
		Vecd getVStar(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
	};

	class DissipativeRiemannSolver
	{
		Fluid& fluid_i_, & fluid_j_;
	public:
		DissipativeRiemannSolver(Fluid& fluid_i, Fluid& fluid_j) :
			fluid_i_(fluid_i), fluid_j_(fluid_j) {};
		Real getPStar(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
		Vecd getVStar(const FluidState& state_i, const FluidState& state_j, const Vecd& direction_to_i);
	};
}
