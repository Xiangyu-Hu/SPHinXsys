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
 * @file 	riemann_solvers_ck.h
 * @brief 	This is the collection of Riemann solvers for weakly compressible fluids.
 * @author	Xiangyu Hu
 */

#ifndef RIEMANN_SOLVER_CK_H
#define RIEMANN_SOLVER_CK_H

#include "data_type.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
template <typename...>
class ImpedanceModel;

template <class FluidI, class FluidJ>
class ImpedanceModel<FluidI, FluidJ>
{
    Real rho0_i_, rho0_j_, rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_;
    Real inv_rho0c0_ave_, rho0c0_geo_ave_, inv_c0_ave_;

  public:
    ImpedanceModel(FluidI &fluid_i, FluidJ &fluid_j);
    Real Rho_i(UnsignedInt) const { return rho0_i_; };
    Real Rho_j(UnsignedInt) const { return rho0_j_; };
    Real Impedance_i(UnsignedInt) const { return rho0c0_i_; };
    Real Impedance_j(UnsignedInt) const { return rho0c0_j_; };
    Real InvImpedanceSum(UnsignedInt, UnsignedInt) const { return inv_rho0c0_sum_; };
    Real InvImpedanceAve(UnsignedInt, UnsignedInt) const { return inv_rho0c0_ave_; };
    Real ImpedanceGeoAve(UnsignedInt, UnsignedInt) const { return rho0c0_geo_ave_; };
    Real InvSoundSpeedAve(UnsignedInt, UnsignedInt) const { return inv_c0_ave_; };
};

class NotUsed;
class Base;

template <typename...>
class RiemannSolver;

template <class FluidI, class FluidJ>
class RiemannSolver<Base, FluidI, FluidJ> : public ImpedanceModel<FluidI, FluidJ>
{
  public:
    typedef FluidI SourceFluid;
    typedef FluidJ TargetFluid;
    RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j);

    template <typename T>
    T AverageP(UnsignedInt i, UnsignedInt j, const T &p_i, const T &p_j) const;
    Vecd AverageV(UnsignedInt i, UnsignedInt j, const Vecd &vel_i, const Vecd &vel_j) const;
};

template <class FluidI, class FluidJ>
class RiemannSolver<NotUsed, FluidI, FluidJ> : public RiemannSolver<Base, FluidI, FluidJ>
{
  public:
    RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
        : RiemannSolver<Base, FluidI, FluidJ>(fluid_i, fluid_j){};
    Real DissipativePJump(UnsignedInt i, UnsignedInt j, const Real &u_jump) const { return 0.0; };
    Real DissipativeUJump(UnsignedInt i, UnsignedInt j, const Real &p_jump) const { return 0.0; };
};
using NoRiemannSolverCK = RiemannSolver<NotUsed, WeaklyCompressibleFluid, WeaklyCompressibleFluid>;

template <class FluidI, class FluidJ, typename LimiterType>
class RiemannSolver<FluidI, FluidJ, LimiterType> : public RiemannSolver<NotUsed, FluidI, FluidJ>
{
  public:
    RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_coeff = 3.0);
    Real DissipativePJump(UnsignedInt i, UnsignedInt j, const Real &u_jump) const;
    Real DissipativeUJump(UnsignedInt i, UnsignedInt j, const Real &p_jump) const;

  protected:
    LimiterType limiter_;
    Real limiter_coeff_;
};
using AcousticRiemannSolverCK = RiemannSolver<WeaklyCompressibleFluid, WeaklyCompressibleFluid, TruncatedLinear>;
using DissipativeRiemannSolverCK = RiemannSolver<WeaklyCompressibleFluid, WeaklyCompressibleFluid, NoLimiter>;
} // namespace SPH
#endif // RIEMANN_SOLVER_CK_H
