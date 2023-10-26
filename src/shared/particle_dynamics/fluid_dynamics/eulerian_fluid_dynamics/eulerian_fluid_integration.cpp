#include "eulerian_fluid_integration.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
FluidStarState EulerianAcousticRiemannSolver::
    getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
    Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
    Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);

    Real p_star = (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ +
                   rhol_cl * rhor_cr * (ul - ur) * SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1))) /
                  (rhol_cl + rhor_cr);
    Real u_star = (rhol_cl * ul + rhor_cr * ur +
                   (state_i.p_ - state_j.p_) * pow(SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1)), 2)) /
                  (rhol_cl + rhor_cr);
    Vecd vel_star = (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_) -
                    e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));

    FluidStarState interface_state(vel_star, p_star);
    interface_state.vel_ = vel_star;
    interface_state.p_ = p_star;

    return interface_state;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
