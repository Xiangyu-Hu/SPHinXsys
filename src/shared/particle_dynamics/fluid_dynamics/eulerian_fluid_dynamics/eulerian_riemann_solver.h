/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file eulerian_riemann_solvers.h
 * @brief This is the collection of Riemann solvers for Eulerian fluid dynamics.
 * @author Zhentong Wang, Chi Zhang and Xiangyu Hu
 */

#ifndef EULERIAN_RIEMANN_SOLVER_H
#define EULERIAN_RIEMANN_SOLVER_H

#include "base_data_type_package.h"
#include "compressible_fluid.h"
#include "riemann_solver.h"
namespace SPH
{
/**
 * @struct CompressibleFluidState
 * @brief  Struct for stored states of Riemann solver in compressible flow.
 */
struct CompressibleFluidState : FluidStateIn
{
    Real &E_;
    CompressibleFluidState(Real &rho, Vecd &vel, Real &p, Real &E)
        : FluidStateIn(rho, vel, p), E_(E) {};
};
struct CompressibleFluidStarState : FluidStateOut
{
    Real E_;
    CompressibleFluidStarState(Real rho, Vecd vel, Real p, Real E)
        : FluidStateOut(rho, vel, p), E_(E) {};
};

/**
 * @struct NoRiemannSolverInCompressibleEulerianMethod
 * @brief  NO RiemannSolver for weakly-compressible flow in Eulerian method for compressible flow.
 */
class NoRiemannSolverInCompressibleEulerianMethod
{
    CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

  public:
    NoRiemannSolverInCompressibleEulerianMethod(CompressibleFluid &fluid_i, CompressibleFluid &fluid_j);
    CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij);
};

/**
 * @struct HLLCRiemannSolver
 * @brief  HLLC Riemann solver.
 */
class HLLCRiemannSolver
{
    CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

  public:
    HLLCRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j, Real limiter_parameter = 0.0);
    CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij);
};
/**
 * @struct HLLCWithLimiterRiemannSolver
 * @brief  HLLC Riemann solver with dissipation limiter.
 */
class HLLCWithLimiterRiemannSolver
{
    CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;
    Real limiter_parameter_;

  public:
    HLLCWithLimiterRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j, Real limiter_parameter = 1.0);
    CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij);
};
} // namespace SPH
#endif // EULERIAN_RIEMANN_SOLVER_H
