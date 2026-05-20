#ifndef RIEMANN_SOLVER_CK_HPP
#define RIEMANN_SOLVER_CK_HPP

#include "riemann_solver_ck.h"

namespace SPH
{
//=================================================================================================//
template <class FluidI, class FluidJ>
RiemannSolver<Base, FluidI, FluidJ>::RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
    : ImpedanceModel<FluidI, FluidJ>(fluid_i, fluid_j) {}
//=================================================================================================//
template <class FluidI, class FluidJ>
template <typename T>
T RiemannSolver<Base, FluidI, FluidJ>::AverageP(
    UnsignedInt i, UnsignedInt j, const T &p_i, const T &p_j) const
{
    return this->InvImpedanceSum(i, j) * (p_i * this->Impedance_j(j) + p_j * this->Impedance_i(i));
}
//=================================================================================================//
template <class FluidI, class FluidJ>
Vecd RiemannSolver<Base, FluidI, FluidJ>::AverageV(
    UnsignedInt i, UnsignedInt j, const Vecd &vel_i, const Vecd &vel_j)
{
    return this->InvImpedanceSum(i, j) * (vel_i * this->Impedance_i(i) + vel_j * this->Impedance_j(j));
}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
RiemannSolver<FluidI, FluidJ, LimiterType>::
    RiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_coeff)
    : RiemannSolver<NotUsed, FluidI, FluidJ>(fluid_i, fluid_j),
      limiter_(this->SoundSpeedAve(0, 0), limiter_coeff) {}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
Real RiemannSolver<FluidI, FluidJ, LimiterType>::DissipativePJump(
    UnsignedInt i, UnsignedInt j, const Real &u_jump)
{
    return this->ImpedanceGeoAve(i, j) *
           u_jump * limiter_(SMAX(u_jump, Real(0))); // the factor 0.5 canceled
}
//=================================================================================================//
template <class FluidI, class FluidJ, typename LimiterType>
Real RiemannSolver<FluidI, FluidJ, LimiterType>::DissipativeUJump(
    UnsignedInt i, UnsignedInt j, const Real &p_jump)
{
    return p_jump * this->InvImpedanceAve(i, j); // the factor 0.5 canceled
}
//=================================================================================================//
} // namespace SPH
#endif // RIEMANN_SOLVER_CK_HPP
