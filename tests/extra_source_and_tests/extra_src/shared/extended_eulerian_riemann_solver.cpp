#include "extended_eulerian_riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class FluidI, class FluidJ>
ExtendedHLLCRiemannSolver::ExtendedHLLCRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_parameter)
    : AcousticRiemannSolver(fluid_i, fluid_j, limiter_parameter), fluid_i_(fluid_i), fluid_j_(fluid_j), limiter_parameter_(limiter_parameter){};
//=================================================================================================//
ExtendedFluidStarState ExtendedHLLCRiemannSolver::
    getExtendedInterfaceState(const ExendedFluidState &state_i, const ExendedFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = ul - fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real s_r = ur + fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) /
                  (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));

    FluidStateIn fluid_state_i(state_i.rho_, state_i.vel_, state_i.p_);
    FluidStateIn fluid_state_j(state_j.rho_, state_j.vel_, state_j.p_);
    FluidStateOut star_state = AcousticRiemannSolver::InterfaceState(fluid_state_i, fluid_state_j, e_ij);

    Real k_star = 0.0;
    Real eps_star = 0.0;

    if (0.0 < s_l)
    {
        k_star = state_i.K_;
        eps_star = state_i.Eps_;

    }
    if (s_l <= 0.0 && 0.0 <= s_star)
    {   
        k_star = state_i.K_ * (s_l - ul) / (s_l - s_star);
        eps_star = state_i.Eps_ * (s_l - ul) / (s_l - s_star);
    }
    if (s_star <= 0.0 && 0.0 <= s_r)
    {
        k_star = state_j.K_ * (s_r - ur) / (s_r - s_star);
        eps_star = state_j.Eps_ * (s_r - ur) / (s_r - s_star);
        
        
    }
    if (s_r < 0.0)
    {
        k_star = state_j.K_;
        eps_star = state_j.Eps_;
    }
    return ExtendedFluidStarState(
        star_state.rho_,
        star_state.vel_,
        star_state.p_,
        k_star,
        eps_star);
}
//=================================================================================================//
template <class FluidI, class FluidJ>
SecondOrderUpwind::SecondOrderUpwind(FluidI &fluid_i, FluidJ &fluid_j, Real limiter_parameter)
    : AcousticRiemannSolver(fluid_i, fluid_j, limiter_parameter),
      fluid_i_(fluid_i), fluid_j_(fluid_j), limiter_parameter_(limiter_parameter){};
//=================================================================================================//
ExtendedFluidStarState SecondOrderUpwind::
    getExtendedInterfaceState(const FluidStateSecondOrderUpwind &state_i, const FluidStateSecondOrderUpwind &state_j, const Vecd &e_ij, const Real &r_ij)
{
    FluidStateIn fluid_state_i(state_i.rho_, state_i.vel_, state_i.p_);
    FluidStateIn fluid_state_j(state_j.rho_, state_j.vel_, state_j.p_);
    FluidStateOut star_state = AcousticRiemannSolver::InterfaceState(fluid_state_i, fluid_state_j, e_ij);

    Real vel_dot_eij = star_state.vel_.dot(e_ij);
    Real phi = 0.0;
    Real r = 0.0;
    Vecd K_grad_limited = Vecd::Zero();
    Vecd Eps_grad_limited = Vecd::Zero();
    Real k_star = 0.0;
    Real eps_star = 0.0;
    if (vel_dot_eij > 0) // Flow from j to i (j is upstream)
    {
        // Upstream: K_j, Downstream: K_i
        // r = (K_grad_[index_j].dot(e_ij)) / (K_grad_[index_i].dot(e_ij));
        // Apply limiter (Minmod)
        // phi = std::max(0.0, std::min(1.0, r));
        K_grad_limited = state_j.K_grad_;
        k_star = state_j.K_ + 0.5 * K_grad_limited.dot(e_ij) * r_ij;
        Eps_grad_limited = state_j.Eps_grad_;
        eps_star = state_j.Eps_ + 0.5 * Eps_grad_limited.dot(e_ij) * r_ij;
    }
    else // Flow from i to j (i is upstream)
    {
        // Upstream: K_i, Downstream: K_j
        // r = (K_grad_[index_i].dot(-e_ij)) / (K_grad_[index_j].dot(-e_ij)); // Ratio-based r for upstream and downstream
        // Apply limiter (e.g., Minmod)
        // phi = std::max(0.0, std::min(1.0, r));
        K_grad_limited = state_i.K_grad_;
        // Compute interface value
        k_star = state_i.K_ + 0.5 * K_grad_limited.dot(-e_ij) * r_ij;
        Eps_grad_limited = state_i.Eps_grad_;
        eps_star = state_i.Eps_ + 0.5 * Eps_grad_limited.dot(-e_ij) * r_ij;
    }
    return ExtendedFluidStarState(
        star_state.rho_,
        star_state.vel_,
        star_state.p_,
        k_star,
        eps_star);
}
} // namespace fluid_dynamics
} // namespace SPH