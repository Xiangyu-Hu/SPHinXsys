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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
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

#ifndef EXTENDED_EULERIAN_RIEMANN_SOLVER_H
#define EXTENDED_EULERIAN_RIEMANN_SOLVER_H

#include "base_data_package.h"
#include "compressible_fluid.h"
#include "eulerian_riemann_solver.h"
namespace SPH
{

struct ExendedFluidState : FluidStateIn
{
    Real &K_;
    Real &Eps_;
    ExendedFluidState(Real &rho, Vecd &vel, Real &p, Real &K, Real &Eps)
        : FluidStateIn(rho, vel, p), K_(K), Eps_(Eps){};
};

struct ExtendedFluidStarState : FluidStateOut
{
    Real K_;
    Real Eps_; 
    ExtendedFluidStarState(Real rho, Vecd vel, Real p, Real K, Real Eps)
        : FluidStateOut(rho, vel, p), K_(K), Eps_(Eps){};
};

class ExtendedHLLCRiemannSolver : NoRiemannSolver
{
    Fluid &fluid_i_, &fluid_j_;
  public:
  template <class FluidI, class FluidJ>
    ExtendedHLLCRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j);
    ExtendedFluidStarState getExtendedInterfaceState(const ExendedFluidState &state_i, const ExendedFluidState &state_j, const Vecd &e_ij);
};
} // namespace SPH
#endif // EXTENDED_EULERIAN_RIEMANN_SOLVER_H