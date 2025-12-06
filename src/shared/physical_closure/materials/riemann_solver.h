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
 * @file 	riemann_solvers.h
 * @brief 	This is the collection of Riemann solvers.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include "base_data_type_package.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
struct FluidStateIn
{
    Vecd &vel_;
    Real &rho_, &p_;
    FluidStateIn(Real &rho, Vecd &vel, Real &p) : vel_(vel), rho_(rho), p_(p) {};
};

struct FluidStateOut
{
    Vecd vel_;
    Real rho_, p_;
    FluidStateOut(Real rho, Vecd vel, Real p) : vel_(vel), rho_(rho), p_(p) {};
};

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

    template <typename T>
    T AverageP(const T &p_i, const T &p_j)
    {
        return (p_i * rho0c0_j_ + p_j * rho0c0_i_) * inv_rho0c0_sum_;
    };

    Vecd AverageV(const Vecd &vel_i, const Vecd &vel_j);
    FluidStateOut InterfaceState(const FluidStateIn &state_i, const FluidStateIn &state_j, const Vecd &e_ij);

  protected:
    Real rho0_i_, rho0_j_;
    Real c0_i_, c0_j_;
    Real rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_;
};

template <typename LimiterType>
class BaseAcousticRiemannSolver : public NoRiemannSolver
{
  public:
    template <class FluidI, class FluidJ>
    BaseAcousticRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_coeff = 3.0)
        : NoRiemannSolver(fluid_i, fluid_j),
          inv_rho0c0_ave_((rho0c0_i_ + rho0c0_j_) / (math::pow(rho0c0_i_, 2) + math::pow(rho0c0_j_, 2))),
          rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
          limiter_(0.5 * (rho0_i_ + rho0_j_) * inv_rho0c0_ave_, limiter_coeff){};
    Real DissipativePJump(const Real &u_jump)
    {
        return rho0c0_geo_ave_ * u_jump * limiter_(SMAX(u_jump, Real(0)));
    };
    Real DissipativeUJump(const Real &p_jump)
    {
        return p_jump * inv_rho0c0_ave_;
    };

    FluidStateOut InterfaceState(const FluidStateIn &state_i, const FluidStateIn &state_j, const Vecd &e_ij)
    {
        FluidStateOut average_state = NoRiemannSolver::InterfaceState(state_i, state_j, e_ij);

        Real ul = -e_ij.dot(state_i.vel_);
        Real ur = -e_ij.dot(state_j.vel_);
        Real u_jump = ul - ur;
        Real limited_mach_number = limiter_(SMAX(u_jump, Real(0)));

        Real p_star = average_state.p_ + 0.5 * rho0c0_geo_ave_ * u_jump * limited_mach_number;
        Real u_dissipative = 0.5 * (state_i.p_ - state_j.p_) * inv_rho0c0_ave_ * limited_mach_number * limited_mach_number;
        Vecd vel_star = average_state.vel_ - e_ij * u_dissipative;

        return FluidStateOut(average_state.rho_, vel_star, p_star);
    };

  protected:
    Real inv_rho0c0_ave_, rho0c0_geo_ave_;
    LimiterType limiter_;
    Real limiter_coeff_;
};
using AcousticRiemannSolver = BaseAcousticRiemannSolver<TruncatedLinear>;
using DissipativeRiemannSolver = BaseAcousticRiemannSolver<NoLimiter>;

} // namespace SPH
#endif // RIEMANN_SOLVER_H
