/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	riemann_solvers.h
 * @brief 	This is the collection of Riemann solvers. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include "base_data_package.h"

namespace SPH
{
	/**
	 * @struct FluidState
	 * @brief  Struct for stored states of Riemann solver in weakly-compressible flow. 
	 */
	struct FluidState
	{
		Vecd &vel_;
		Real &rho_, &p_;
		FluidState(Real &rho, Vecd &vel, Real &p)
			: vel_(vel), rho_(rho), p_(p){};
	};

	struct FluidStarState
	{
		Vecd vel_;
		Real p_;
		FluidStarState(Vecd vel, Real p)
			: vel_(vel), p_(p){};
	};

	/**
	 * @struct CompressibleFluidState
	 * @brief  Struct for stored states of Riemann solver in compressible flow. 
	 */

	struct CompressibleFluidState : FluidState
	{
		Real &E_;
		CompressibleFluidState(Real &rho, Vecd &vel, Real &p, Real &E)
			: FluidState(rho, vel, p), E_(E){};
	};

	struct CompressibleFluidStarState : FluidStarState
	{
		Real rho_;
		Real E_;
		CompressibleFluidStarState(Real rho, Vecd vel, Real p, Real E)
			: FluidStarState(vel, p), rho_(rho), E_(E){};
	};

	class Fluid;
	class CompressibleFluid;

	/**
	 * @struct NoRiemannSolver
	 * @brief  Central difference scheme without Riemann flux. 
	 */
	class NoRiemannSolver
	{
	public:
		template <class FluidI, class FluidJ>
		NoRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
			: rho0_i_(fluid_i.ReferenceDensity()), rho0_j_(fluid_j.ReferenceDensity()),
			  c0_i_(fluid_i.ReferenceSoundSpeed()), c0_j_(fluid_j.ReferenceSoundSpeed()),
			  rho0c0_i_(rho0_i_ * c0_i_), rho0c0_j_(rho0_j_ * c0_j_),
			  inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)){};
		Real DissipativePJump(const Real &u_jump);
		Real DissipativeUJump(const Real &p_jump);
		Real AverageP(const Real &p_i, const Real &p_j);
		Vecd AverageV(const Vecd &vel_i, const Vecd &vel_j);

	protected:
		Real rho0_i_, rho0_j_;
		Real c0_i_, c0_j_;
		Real rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_;
	};

	class AcousticRiemannSolver : public NoRiemannSolver
	{
	public:
		template <class FluidI, class FluidJ>
		AcousticRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
			: NoRiemannSolver(fluid_i, fluid_j),
			  inv_rho0c0_ave_(2.0 * inv_rho0c0_sum_),
			  rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
			  inv_c_ave_(0.5 * (rho0_i_ + rho0_j_) * inv_rho0c0_ave_){};
		Real DissipativePJump(const Real &u_jump);
		Real DissipativeUJump(const Real &p_jump);

	protected:
		Real inv_rho0c0_ave_, rho0c0_geo_ave_;
		Real inv_c_ave_;
	};

	class DissipativeRiemannSolver : public AcousticRiemannSolver
	{
	public:
		template <class FluidI, class FluidJ>
		DissipativeRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
			: AcousticRiemannSolver(fluid_i, fluid_j){};
		Real DissipativePJump(const Real &u_jump);
	};

	/**
	 * @struct HLLCRiemannSolverInWeaklyCompressibleFluid
	 * @brief  HLLC Riemann for weakly-compressible flow. 
	 */
	class HLLCRiemannSolverInWeaklyCompressibleFluid
	{
		Fluid &fluid_i_, &fluid_j_;

	public:
		HLLCRiemannSolverInWeaklyCompressibleFluid(Fluid &compressible_fluid_i, Fluid &compressible_fluid_j)
			: fluid_i_(compressible_fluid_i), fluid_j_(compressible_fluid_j){};
		FluidState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
	};

	/**
	 * @struct HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid
	 * @brief  HLLC Riemann with dissipation limiter for weakly-compressible flow. 
	 */
	class HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid
	{
		Fluid &fluid_i_, &fluid_j_;

	public:
		HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid(Fluid &compressible_fluid_i, Fluid &compressible_fluid_j) : fluid_i_(compressible_fluid_i), fluid_j_(compressible_fluid_j){};
		FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &direction_to_i);
	};

	/**
	 * @struct HLLCRiemannSolver
	 * @brief  HLLC Riemann solver. 
	 */
	class HLLCRiemannSolver
	{
		CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

	public:
		HLLCRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j)
			: compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
		CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
	};

	/**
	 * @struct HLLCRiemannSolver
	 * @brief  HLLC Riemann solver with dissipation limiter. 
	 */
	class HLLCWithLimiterRiemannSolver
	{
		CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

	public:
		HLLCWithLimiterRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j)
			: compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
		CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &direction_to_i);
	};
}

#endif // RIEMANN_SOLVER_H
