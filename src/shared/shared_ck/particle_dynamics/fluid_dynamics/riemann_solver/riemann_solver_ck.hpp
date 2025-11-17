#ifndef RIEMANN_SOLVER_CK_HPP
#define RIEMANN_SOLVER_CK_HPP

#include "riemann_solver_ck.h"

namespace SPH
{
//=================================================================================================//
template <class FluidI, class FluidJ>
RiemannSolver<Base, FluidI, FluidJ>::RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
    : rho0_i_(fluid_i.ReferenceDensity()), rho0_j_(fluid_j.ReferenceDensity()),
      c0_i_(fluid_i.ReferenceSoundSpeed()), c0_j_(fluid_j.ReferenceSoundSpeed()),
      rho0c0_i_(rho0_i_ * c0_i_), rho0c0_j_(rho0_j_ * c0_j_),
      inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)) {}
//=================================================================================================//
template <class FluidI, class FluidJ>
template <typename T>
T RiemannSolver<Base, FluidI, FluidJ>::AverageP(const T &p_i, const T &p_j) const
{
    return (p_i * rho0c0_j_ + p_j * rho0c0_i_) * inv_rho0c0_sum_;
}
//=================================================================================================//
template <class FluidI, class FluidJ>
Vecd RiemannSolver<Base, FluidI, FluidJ>::AverageV(const Vecd &vel_i, const Vecd &vel_j)
{
    return (vel_i * rho0c0_i_ + vel_j * rho0c0_j_) * inv_rho0c0_sum_;
}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
RiemannSolver<FluidI, FluidJ, LimiterType>::
    RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_coeff)
    : RiemannSolver<NotUsed, FluidI, FluidJ>(fluid_i, fluid_j),
      inv_rho0c0_ave_((this->rho0c0_i_ + this->rho0c0_j_) /
                      (math::pow(this->rho0c0_i_, 2) + math::pow(this->rho0c0_j_, 2))),
      rho0c0_geo_ave_(2.0 * this->rho0c0_i_ * this->rho0c0_j_ * this->inv_rho0c0_sum_),
      limiter_(0.5 * (this->rho0_i_ + this->rho0_j_) * inv_rho0c0_ave_, limiter_coeff) {}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
Real RiemannSolver<FluidI, FluidJ, LimiterType>::DissipativePJump(const Real &u_jump)
{
    return rho0c0_geo_ave_ * u_jump * limiter_(SMAX(u_jump, Real(0))); // the factor 0.5 canceled
}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
Real RiemannSolver<FluidI, FluidJ, LimiterType>::DissipativeUJump(const Real &p_jump)
{
    return p_jump * inv_rho0c0_ave_; // the factor 0.5 canceled with 2.0 in kernel approximation
}
//=================================================================================================//
} // namespace SPH
#endif // RIEMANN_SOLVER_CK_HPP
