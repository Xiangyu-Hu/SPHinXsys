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
 * @file 	riemann_solvers.h
 * @brief 	This is the collection of Riemann solvers.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include "base_data_package.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
struct FluidStateIn
{
    Vecd &vel_;
    Real &rho_, &p_;
    FluidStateIn(Real &rho, Vecd &vel, Real &p) : vel_(vel), rho_(rho), p_(p){};
};

struct FluidStateOut
{
    Vecd vel_;
    Real rho_, p_;
    FluidStateOut(Real rho, Vecd vel, Real p) : vel_(vel), rho_(rho), p_(p){};
};

class NoRiemannSolverKernel
{
  public:
    NoRiemannSolverKernel(DeviceReal rho0c0_i, DeviceReal rho0c0_j, DeviceReal inv_rho0c0_sum)
        : rho0c0_i_(rho0c0_i), rho0c0_j_(rho0c0_j), inv_rho0c0_sum_(inv_rho0c0_sum) {}

    DeviceReal DissipativePJump(const DeviceReal &) const { return DeviceReal(0.0); }
    DeviceReal DissipativeUJump(const DeviceReal &) const { return DeviceReal(0.0); }

    DeviceReal AverageP(const DeviceReal &p_i, const DeviceReal &p_j)
    {
        return (p_i * rho0c0_j_ + p_j * rho0c0_i_) * inv_rho0c0_sum_;
    };
  protected:
    DeviceReal rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_;
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
          inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)),
          device_kernel(rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_) {};
    Real DissipativePJump(const Real &u_jump) const { return 0.0; }
    Real DissipativeUJump(const Real &p_jump) const { return 0.0; }

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

  public:
    execution::DeviceImplementation<NoRiemannSolverKernel> device_kernel;
};

class AcousticRiemannSolverKernel : public NoRiemannSolverKernel
{
  public:
    template<class... Args>
    AcousticRiemannSolverKernel(DeviceReal inv_rho0c0_ave, DeviceReal rho0c0_geo_ave,
                                DeviceReal inv_c_ave, DeviceReal limiter_coeff, Args&& ...baseArgs)
        : NoRiemannSolverKernel(std::forward<Args>(baseArgs)...), inv_rho0c0_ave_(inv_rho0c0_ave),
          rho0c0_geo_ave_(rho0c0_geo_ave), inv_c_ave_(inv_c_ave), limiter_coeff_(limiter_coeff) {}

    DeviceReal DissipativePJump(const DeviceReal &u_jump) const {
        return rho0c0_geo_ave_ * u_jump * sycl::min(limiter_coeff_ * sycl::max(u_jump * inv_c_ave_, DeviceReal(0)), DeviceReal(1));
    }

    DeviceReal DissipativeUJump(const DeviceReal &p_jump) const  {
        return p_jump * inv_rho0c0_ave_;
    }

  protected:
    DeviceReal inv_rho0c0_ave_, rho0c0_geo_ave_;
    DeviceReal inv_c_ave_;
    DeviceReal limiter_coeff_;
};

class AcousticRiemannSolver : public NoRiemannSolver
{
  public:
    template <class FluidI, class FluidJ>
    AcousticRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, const Real limiter_coeff = 3.0)
        : NoRiemannSolver(fluid_i, fluid_j),
          inv_rho0c0_ave_(2.0 * inv_rho0c0_sum_),
          rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
          inv_c_ave_(0.5 * (rho0_i_ + rho0_j_) * inv_rho0c0_ave_),
          limiter_coeff_(limiter_coeff),
          device_kernel(inv_rho0c0_ave_, rho0c0_geo_ave_, inv_c_ave_, limiter_coeff_,
                        rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_) {};
    Real DissipativePJump(const Real &u_jump) const {
        return rho0c0_geo_ave_ * u_jump * SMIN(limiter_coeff_ * SMAX(u_jump * inv_c_ave_, Real(0)), Real(1));
    }
    Real DissipativeUJump(const Real &p_jump) const  {
        return p_jump * inv_rho0c0_ave_;
    }
    FluidStateOut InterfaceState(const FluidStateIn &state_i, const FluidStateIn &state_j, const Vecd &e_ij);

  protected:
    Real inv_rho0c0_ave_, rho0c0_geo_ave_;
    Real inv_c_ave_;
    Real limiter_coeff_;

  public:
    execution::DeviceImplementation<AcousticRiemannSolverKernel> device_kernel;
};

class DissipativeRiemannSolver : public AcousticRiemannSolver
{
  public:
    template <class FluidI, class FluidJ>
    DissipativeRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
        : AcousticRiemannSolver(fluid_i, fluid_j){};
    Real DissipativePJump(const Real &u_jump) const {
        return rho0c0_geo_ave_ * u_jump;
    }
};
} // namespace SPH

#endif // RIEMANN_SOLVER_H
