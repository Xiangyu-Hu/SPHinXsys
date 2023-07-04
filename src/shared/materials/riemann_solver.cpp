#include "riemann_solver.h"
#include "base_material.h"

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
Real NoRiemannSolver::AverageP(const Real &p_i, const Real &p_j)
{
    return (p_i * rho0c0_j_ + p_j * rho0c0_i_) * inv_rho0c0_sum_;
}
//=================================================================================================//
Vecd NoRiemannSolver::AverageV(const Vecd &vel_i, const Vecd &vel_j)
{
    return (vel_i * rho0c0_i_ + vel_j * rho0c0_j_) * inv_rho0c0_sum_;
}
//=================================================================================================//
Real AcousticRiemannSolver::DissipativePJump(const Real &u_jump)
{
    return rho0c0_geo_ave_ * u_jump * SMIN(Real(3) * SMAX(u_jump * inv_c_ave_, Real(0)), Real(1));
}
//=================================================================================================//
Real AcousticRiemannSolver::DissipativeUJump(const Real &p_jump)
{
    return p_jump * inv_rho0c0_ave_;
}
//=================================================================================================//
Real DissipativeRiemannSolver::DissipativePJump(const Real &u_jump)
{
    return rho0c0_geo_ave_ * u_jump;
}
//=================================================================================================//
} // namespace SPH
