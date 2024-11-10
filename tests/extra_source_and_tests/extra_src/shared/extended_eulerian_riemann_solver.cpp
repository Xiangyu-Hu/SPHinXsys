#include "extended_eulerian_riemann_solver.h"

namespace SPH
{
//=================================================================================================//
template <class FluidI, class FluidJ>
ExtendedHLLCRiemannSolver::ExtendedHLLCRiemannSolver(FluidI &fluid_i, FluidJ &fluid_j)
    : NoRiemannSolver(fluid_i, fluid_j), fluid_i_(fluid_i), fluid_j_(fluid_j)
    {};
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
Real p_star = 0.0;
Vecd v_star = Vecd::Zero();
Real rho_star = 0.0;
Real k_star = 0.0;
Real eps_star = 0.0;

if (0.0 < s_l)
{
    p_star = state_i.p_;
    v_star = state_i.vel_;
    rho_star = state_i.rho_;
    k_star = state_i.K_; // Use the left state's TKE
    eps_star = state_i.Eps_;
    
}
if (s_l <= 0.0 && 0.0 <= s_star)
{
    p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
    v_star = state_i.vel_ - e_ij * (s_star - ul);
    rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
    k_star = state_i.K_ * (s_l - ul) / (s_l - s_star);
    eps_star = state_i.Eps_ * (s_l - ul) / (s_l - s_star);
   
}
if (s_star <= 0.0 && 0.0 <= s_r)
{
    p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
    v_star = state_j.vel_ - e_ij * (s_star - ur);
    rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
    k_star = state_j.K_ * (s_r - ur) / (s_r - s_star);
    eps_star = state_j.Eps_ * (s_r - ur) / (s_r - s_star);

}
if (s_r < 0.0)
{
    p_star = state_j.p_;
    v_star = state_j.vel_;
    rho_star = state_j.rho_;
    k_star = state_j.K_;
    eps_star = state_j.Eps_;
}
     return ExtendedFluidStarState(
        rho_star,
        v_star,
        p_star,
        k_star,
        eps_star
    );
}
} // namespace SPH
//=================================================================================================//