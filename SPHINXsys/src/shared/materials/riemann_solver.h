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

#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include "base_data_package.h"

namespace SPH
{
	struct FluidState
	{
		Vecd &vel_;
		Real &rho_, &p_;
		FluidState(Real &rho, Vecd &vel, Real &p)
			: vel_(vel), rho_(rho), p_(p){};
	};

	struct CompressibleFluidState : FluidState
	{
		Real &E_;
		CompressibleFluidState(Real &rho, Vecd &vel, Real &p, Real &E)
			: FluidState(rho, vel, p), E_(E){};
	};

	class Fluid;
	class CompressibleFluid;

	class NoRiemannSolver
	{
		Fluid &fluid_l_, &fluid_r_;

	public:
		NoRiemannSolver(Fluid &fluid_i, Fluid &fluid_j) : fluid_l_(fluid_i), fluid_r_(fluid_j){};
		Real getPStar(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
		Vecd getVStar(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
	};

	class BaseAcousticRiemannSolver
	{
	protected:
		Fluid &fluid_i_, &fluid_j_;

	public:
		BaseAcousticRiemannSolver(Fluid &fluid_i, Fluid &fluid_j) : fluid_i_(fluid_i), fluid_j_(fluid_j){};
		inline void prepareSolver(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i,
								  Real &ul, Real &ur, Real &rhol_cl, Real &rhor_cr);
	};
	class AcousticRiemannSolver : public BaseAcousticRiemannSolver
	{
	public:
		AcousticRiemannSolver(Fluid &fluid_i, Fluid &fluid_j) : BaseAcousticRiemannSolver(fluid_i, fluid_j){};
		Real getPStar(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
		Vecd getVStar(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
	};

	class DissipativeRiemannSolver : public BaseAcousticRiemannSolver
	{
	public:
		DissipativeRiemannSolver(Fluid &fluid_i, Fluid &fluid_j) : BaseAcousticRiemannSolver(fluid_i, fluid_j){};
		Real getPStar(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
		Vecd getVStar(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
	};

	class HLLCRiemannSolver
	{
		CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

	public:
		HLLCRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j)
			: compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
		Real getPStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
		Vecd getVStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
		Real getRhoStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
		Real getEStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
	};

	class HLLCWithLimiterRiemannSolver
	{
		CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

	public:
		HLLCWithLimiterRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j)
			: compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
		Real getPStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
		Vecd getVStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
		Real getRhoStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
		Real getEStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
	};
}

#endif //RIEMANN_SOLVER_H
