#include "riemann_solver.h"
#include "base_material.h"

namespace SPH
{
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
} // namespace SPH
