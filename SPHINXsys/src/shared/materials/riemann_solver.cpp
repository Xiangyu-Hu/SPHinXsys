/**
 * @file 	riemann_solver.cpp
 * @author	Xiangyu Hu
 */

#include "riemann_solver.h"

#include "base_material.h"

namespace SPH {
	//=================================================================================================//
	Real NoRiemannSolver::
		getPStar(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
	{
		return (state_i.p_ * state_j.rho_ + state_j.p_ * state_i.rho_) 
			/ (state_i.rho_ + state_j.rho_);
	}
	//=================================================================================================//
	Vecd NoRiemannSolver::
		getVStar(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
	{
		return (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) 
			/ (state_i.rho_ + state_j.rho_);
	}
	//=================================================================================================//
	Real AcousticRiemannSolver::
		getPStar(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);

		Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
		Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
		Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);

		return (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ + rhol_cl * rhor_cr * (ul - ur)
			* SMIN(3.0 * SMAX((ul - ur) / clr, 0.0), 1.0)) / (rhol_cl + rhor_cr);
	}
	//=================================================================================================//
	Vecd AcousticRiemannSolver::
		getVStar(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
		Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
		Real u_star = (rhol_cl * ul + rhor_cr * ur + state_i.p_ - state_j.p_) / (rhol_cl + rhor_cr);

		return (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_)
			- e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));
	}
	//=================================================================================================//
	Real DissipativeRiemannSolver::
		getPStar(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);

		Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
		Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
		Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);

		return (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ + rhol_cl * rhor_cr * (ul - ur)) / (rhol_cl + rhor_cr);
	}
	//=================================================================================================//
	Vecd DissipativeRiemannSolver::
		getVStar(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
		Real rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
		Real u_star = (rhol_cl * ul + rhor_cr * ur + state_i.p_ - state_j.p_) / (rhol_cl + rhor_cr);

		return (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_)
			- e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));
	}
	//=================================================================================================//
}
