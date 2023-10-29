#include "eulerian_riemann_solver.h"

namespace SPH
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
NoRiemannSolverInCompressibleEulerianMethod::
    NoRiemannSolverInCompressibleEulerianMethod(CompressibleFluid &compressible_fluid_i,
                                                CompressibleFluid &compressible_fluid_j)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
//=================================================================================================//
CompressibleFluidStarState NoRiemannSolverInCompressibleEulerianMethod::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real p_star = 0.5 * (state_i.p_ + state_j.p_);
    Vecd v_star = 0.5 * (state_i.vel_ + state_j.vel_);
    Real rho_star = 0.5 * (state_i.rho_ + state_j.rho_);
    Real energy_star = 0.5 * (state_i.E_ + state_j.E_);

    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
}
//=================================================================================================//
HLLCRiemannSolver::HLLCRiemannSolver(CompressibleFluid &compressible_fluid_i,
                                     CompressibleFluid &compressible_fluid_j, Real limiter_parameter)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
//=================================================================================================//
CompressibleFluidStarState HLLCRiemannSolver::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) /
                  (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
    Real p_star = 0.0;
    Vecd v_star = Vecd::Zero();
    Real rho_star = 0.0;
    Real energy_star = 0.0;
    if (0.0 < s_l)
    {
        p_star = state_i.p_;
        v_star = state_i.vel_;
        rho_star = state_i.rho_;
        energy_star = state_i.E_;
    }
    if (s_l <= 0.0 && 0.0 <= s_star)
    {
        p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
        v_star = state_i.vel_ - e_ij * (s_star - ul);
        rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
        energy_star = state_i.rho_ * (s_l - ul) / (s_l - s_star) *
                      (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)));
    }
    if (s_star <= 0.0 && 0.0 <= s_r)
    {
        p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
        v_star = state_j.vel_ - e_ij * (s_star - ur);
        rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
        energy_star = state_j.rho_ * (s_r - ur) / (s_r - s_star) *
                      (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
    }
    if (s_r < 0.0)
    {
        p_star = state_j.p_;
        v_star = state_j.vel_;
        rho_star = state_j.rho_;
        energy_star = state_j.E_;
    }
    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
}
//=================================================================================================//
HLLCWithLimiterRiemannSolver::
    HLLCWithLimiterRiemannSolver(CompressibleFluid &compressible_fluid_i,
                                 CompressibleFluid &compressible_fluid_j, Real limiter_parameter)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j),
      limiter_parameter_(limiter_parameter){};
//=================================================================================================//
CompressibleFluidStarState HLLCWithLimiterRiemannSolver::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real rhol_cl = compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
    Real rhor_cr = compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
    Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);
    Real s_star = (state_j.p_ - state_i.p_) * pow(SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1)), 2) /
                      (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur)) +
                  (state_i.rho_ * (s_l - ul) * ul - state_j.rho_ * (s_r - ur) * ur) /
                      (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur));
    Real p_star = 0.0;
    Vecd v_star = Vecd::Zero();
    Real rho_star = 0.0;
    Real energy_star = 0.0;
    if (0.0 < s_l)
    {
        p_star = state_i.p_;
        v_star = state_i.vel_;
        rho_star = state_i.rho_;
        energy_star = state_i.E_;
    }
    if (s_l <= 0.0 && 0.0 <= s_star)
    {
        p_star = 0.5 * (state_i.p_ + state_j.p_) +
                 0.5 * (state_i.rho_ * (s_l - ul) * (s_star - ul) + state_j.rho_ * (s_r - ur) * (s_star - ur)) *
                     SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1));
        v_star = state_i.vel_ - e_ij * (s_star - ul);
        rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
        energy_star = ((s_l - ul) * state_i.E_ - state_i.p_ * ul + p_star * s_star) / (s_l - s_star);
    }
    if (s_star <= 0.0 && 0.0 <= s_r)
    {
        p_star = 0.5 * (state_i.p_ + state_j.p_) +
                 0.5 * (state_i.rho_ * (s_l - ul) * (s_star - ul) + state_j.rho_ * (s_r - ur) * (s_star - ur)) *
                     SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1));
        v_star = state_j.vel_ - e_ij * (s_star - ur);
        rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
        energy_star = ((s_r - ur) * state_j.E_ - state_j.p_ * ur + p_star * s_star) / (s_r - s_star);
    }
    if (s_r < 0.0)
    {
        p_star = state_j.p_;
        v_star = state_j.vel_;
        rho_star = state_j.rho_;
        energy_star = state_j.E_;
    }
    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
}
//=================================================================================================//
} // namespace SPH
