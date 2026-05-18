#include "riemann_solver.h"

namespace SPH
{
//=================================================================================================//
Real NoRiemannSolver::DissipativePJump(const Real &u_jump)
{
    return 0.0;
}
//=================================================================================================//
Real NoRiemannSolver::DissipativeUJump(const Real &p_jump)
{
    return 0.0;
}
//=================================================================================================//
Vecd NoRiemannSolver::AverageV(const Vecd &vel_i, const Vecd &vel_j)
{
    return (vel_i * rho0c0_i_ + vel_j * rho0c0_j_) * inv_rho0c0_sum_;
}
//=================================================================================================//
FluidStateOut NoRiemannSolver::
    InterfaceState(const FluidStateIn &state_i, const FluidStateIn &state_j, const Vecd &e_ij)
{
    Real rho_star = 0.5 * (state_i.rho_ + state_j.rho_);
    Real p_star = 0.5 * (state_i.p_ + state_j.p_);
    Vecd v_star = 0.5 * (state_i.vel_ + state_j.vel_);

    return FluidStateOut(rho_star, v_star, p_star);
}
//=================================================================================================//
} // namespace SPH
