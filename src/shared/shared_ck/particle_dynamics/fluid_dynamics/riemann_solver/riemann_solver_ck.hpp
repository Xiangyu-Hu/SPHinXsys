#ifndef RIEMANN_SOLVER_CK_HPP
#define RIEMANN_SOLVER_CK_HPP

#include "riemann_solver_ck.h"

namespace SPH
{
//=================================================================================================//
template <class FluidI, class FluidJ>
template <class ExecutionPolicy, class EncloserType>
RiemannSolver<Base, FluidI, FluidJ>::ComputingKernel::ComputingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : ImpedanceModel<FluidI, FluidJ>(ex_policy, encloser.fluid_i_, encloser.fluid_j_) {}
//=================================================================================================//
template <class FluidI, class FluidJ>
template <typename T>
T RiemannSolver<Base, FluidI, FluidJ>::ComputingKernel::AverageP(
    UnsignedInt i, UnsignedInt j, const T &p_i, const T &p_j) const
{
    return this->InvImpedanceSum(i, j) * (p_i * this->Impedance_j(j) + p_j * this->Impedance_i(i));
}
//=================================================================================================//
template <class FluidI, class FluidJ>
Vecd RiemannSolver<Base, FluidI, FluidJ>::ComputingKernel::AverageV(
    UnsignedInt i, UnsignedInt j, const Vecd &vel_i, const Vecd &vel_j) const
{
    return (vel_i * this->Impedance_i(i) + vel_j * this->Impedance_j(j)) * this->InvImpedanceSum(i, j);
}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
RiemannSolver<FluidI, FluidJ, LimiterType>::
    RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_coeff)
    : BaseRiemannSolver(fluid_i, fluid_j), limiter_coeff_(limiter_coeff) {}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
template <class ExecutionPolicy, class EncloserType>
RiemannSolver<FluidI, FluidJ, LimiterType>::ComputingKernel::ComputingKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseRiemannSolver::ComputingKernel(ex_policy, encloser),
      limiter_(encloser.limiter_coeff_) {}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
Real RiemannSolver<FluidI, FluidJ, LimiterType>::ComputingKernel::DissipativePJump(
    UnsignedInt i, UnsignedInt j, const Real &u_jump) const
{
    return this->ImpedanceGeoAve(i, j) * u_jump *
           limiter_(this->InvSoundSpeedAve(i, j) * SMAX(u_jump, Real(0))); // the factor 0.5 canceled
}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
Real RiemannSolver<FluidI, FluidJ, LimiterType>::ComputingKernel::DissipativeUJump(
    UnsignedInt i, UnsignedInt j, const Real &p_jump) const
{
    return p_jump * this->InvImpedanceAve(i, j); // the factor 0.5 canceled
}
//=================================================================================================//
template <class ExecutionPolicy>
ImpedanceModel<WeaklyCompressibleFluid, WeaklyCompressibleFluid>::ImpedanceModel(
    const ExecutionPolicy &ex_policy,
    const WeaklyCompressibleFluid &fluid_i, const WeaklyCompressibleFluid &fluid_j)
    : rho0_i_(fluid_i.ReferenceDensity()), rho0_j_(fluid_j.ReferenceDensity()),
      rho0c0_i_(rho0_i_ * fluid_i.ReferenceSoundSpeed()),
      rho0c0_j_(rho0_j_ * fluid_j.ReferenceSoundSpeed()),
      inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)),
      inv_rho0c0_ave_((rho0c0_i_ + rho0c0_j_) / (math::pow(rho0c0_i_, 2) + math::pow(rho0c0_j_, 2))),
      rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
      inv_c0_ave_(0.5 * (rho0_i_ + rho0_j_) * inv_rho0c0_ave_) {}
//=================================================================================================//
template <class ExecutionPolicy>
ImpedanceModel<WeaklyCompressibleMixture, WeaklyCompressibleMixture>::ImpedanceModel(
    const ExecutionPolicy &ex_policy,
    const WeaklyCompressibleMixture &fluid_i, const WeaklyCompressibleMixture &fluid_j)
    : rho0_i_(fluid_i.dvReferenceDensity()->DelegatedDataView(ex_policy)),
      rho0_j_(fluid_j.dvReferenceDensity()->DelegatedDataView(ex_policy)),
      c0_i_(fluid_i.ReferenceSoundSpeed()), c0_j_(fluid_j.ReferenceSoundSpeed()) {}
//=================================================================================================//
} // namespace SPH
#endif // RIEMANN_SOLVER_CK_HPP
