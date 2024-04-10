#include "riemann_solver.h"

namespace SPH
{
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
FluidStateOut AcousticRiemannSolver::
    InterfaceState(const FluidStateIn &state_i, const FluidStateIn &state_j, const Vecd &e_ij)
{
    FluidStateOut average_state = NoRiemannSolver::InterfaceState(state_i, state_j, e_ij);

    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real u_jump = ul - ur;
    Real limited_mach_number = SMIN(limiter_coeff_ * SMAX(u_jump * inv_c_ave_, Real(0)), Real(1));

    Real p_star = average_state.p_ + 0.5 * rho0c0_geo_ave_ * u_jump * limited_mach_number;
    Real u_dissipative = 0.5 * (state_i.p_ - state_j.p_) * inv_rho0c0_ave_ * limited_mach_number * limited_mach_number;
    Vecd vel_star = average_state.vel_ - e_ij * u_dissipative;

    return FluidStateOut(average_state.rho_, vel_star, p_star);
}
//=================================================================================================//
} // namespace SPH
