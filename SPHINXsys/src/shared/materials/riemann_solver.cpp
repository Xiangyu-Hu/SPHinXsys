/**
 * @file 	riemann_solver.cpp
 * @author	Xiangyu Hu
 */

#include "riemann_solver.h"

#include "base_material.h"

#include "compressible_fluid.h"

namespace SPH
{
	//=================================================================================================//
	Real NoRiemannSolver::
		getPStar(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		return (state_i.p_ * state_j.rho_ + state_j.p_ * state_i.rho_) / (state_i.rho_ + state_j.rho_);
	}
	//=================================================================================================//
	Vecd NoRiemannSolver::
		getVStar(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		return (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_);
	}
	//=================================================================================================//
	void BaseAcousticRiemannSolver::
		prepareSolver(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij,
					  Real &ul, Real &ur, Real &rhol_cl, Real &rhor_cr)
	{
		ul = dot(-e_ij, state_i.vel_);
		ur = dot(-e_ij, state_j.vel_);
		rhol_cl = fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
		rhor_cr = fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
	}
	//=================================================================================================//
	Real AcousticRiemannSolver::
		getPStar(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		Real ul, ur, rhol_cl, rhor_cr;
		prepareSolver(state_i, state_j, e_ij, ul, ur, rhol_cl, rhor_cr);

		Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);
		return (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ +
				rhol_cl * rhor_cr * (ul - ur) * SMIN(3.0 * SMAX((ul - ur) / clr, 0.0), 1.0)) /
			   (rhol_cl + rhor_cr);
	}
	//=================================================================================================//
	Vecd AcousticRiemannSolver::
		getVStar(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		Real ul, ur, rhol_cl, rhor_cr;
		prepareSolver(state_i, state_j, e_ij, ul, ur, rhol_cl, rhor_cr);

		Real u_star = (rhol_cl * ul + rhor_cr * ur + state_i.p_ - state_j.p_) / (rhol_cl + rhor_cr);
		return (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_) -
			   e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));
	}
	//=================================================================================================//
	Real DissipativeRiemannSolver::
		getPStar(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		Real ul, ur, rhol_cl, rhor_cr;
		prepareSolver(state_i, state_j, e_ij, ul, ur, rhol_cl, rhor_cr);

		return (rhol_cl * state_j.p_ + rhor_cr * state_i.p_ + rhol_cl * rhor_cr * (ul - ur)) / (rhol_cl + rhor_cr);
	}
	//=================================================================================================//
	Vecd DissipativeRiemannSolver::
		getVStar(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		Real ul, ur, rhol_cl, rhor_cr;
		prepareSolver(state_i, state_j, e_ij, ul, ur, rhol_cl, rhor_cr);

		Real u_star = (rhol_cl * ul + rhor_cr * ur + state_i.p_ - state_j.p_) / (rhol_cl + rhor_cr);
		return (state_i.vel_ * state_i.rho_ + state_j.vel_ * state_j.rho_) / (state_i.rho_ + state_j.rho_) -
			   e_ij * (u_star - (ul * state_i.rho_ + ur * state_j.rho_) / (state_i.rho_ + state_j.rho_));
	}
	//=================================================================================================//
	Real HLLCRiemannSolver::
		getPStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Real p_star = 0.0;
		if (0.0 < s_l)
		{
			p_star = state_i.p_;
		}
		if (s_l <= 0.0 && s_r >= 0.0)
		{
			p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
		}
		if (s_r < 0.0)
		{
			p_star = state_j.p_;
		}

		return p_star;
	}
	//=================================================================================================//
	Vecd HLLCRiemannSolver::
		getVStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Vecd v_star(0);
		if (0.0 < s_l)
		{
			v_star = state_i.vel_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			v_star = state_i.vel_ - e_ij * (s_star - ul);
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			v_star = state_j.vel_ - e_ij * (s_star - ur);
		}
		if (s_r < 0.0)
		{
			v_star = state_j.vel_;
		}

		return v_star;
	}
	//=================================================================================================//
	Real HLLCRiemannSolver::
		getRhoStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Real rho_star = 0.0;
		if (0.0 < s_l)
		{
			rho_star = state_i.rho_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
		}
		if (s_r < 0.0)
		{
			rho_star = state_j.rho_;
		}

		return rho_star;
	}
	//=================================================================================================//
	Real HLLCRiemannSolver::
		getEStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Real energy_star = 0.0;
		if (0.0 < s_l)
		{
			energy_star = state_i.E_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			energy_star = state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)));
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			energy_star = state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
		}
		if (s_r < 0.0)
		{
			energy_star = state_j.E_;
		}

		return energy_star;
	}
	//=================================================================================================//
	Real HLLCWithLimiterRiemannSolver::
		getPStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Real p_star = 0.0;
		if (0.0 < s_l)
		{
			p_star = state_i.p_;
		}
		if (s_l <= 0.0 && s_r >= 0.0)
		{
			Real rho_ave = 2 * state_i.rho_ * state_j.rho_ / (state_i.rho_ + state_j.rho_);
			Real rho_cl = state_i.rho_ * compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
			Real rho_cr = state_j.rho_ * compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
			Real rho_clr = (rho_cl * state_i.rho_ + rho_cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);
			p_star = 0.5 * (state_i.p_ + state_j.p_) + 0.5 * (SMIN(3.0 * SMAX(rho_ave * (ul - ur), 0.0), rho_clr) * (ul - ur) + s_star * (rho_cr - rho_cl));
		}
		if (s_r < 0.0)
		{
			p_star = state_j.p_;
		}

		return p_star;
	}
	//=================================================================================================//
	Vecd HLLCWithLimiterRiemannSolver::
		getVStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Vecd v_star(0);
		if (0.0 < s_l)
		{
			v_star = state_i.vel_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			v_star = state_i.vel_ - e_ij * (s_star - ul);
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			v_star = state_j.vel_ - e_ij * (s_star - ur);
		}
		if (s_r < 0.0)
		{
			v_star = state_j.vel_;
		}

		return v_star;
	}
	//=================================================================================================//
	Real HLLCWithLimiterRiemannSolver::
		getRhoStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Real rho_star = 0.0;
		if (0.0 < s_l)
		{
			rho_star = state_i.rho_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
		}
		if (s_r < 0.0)
		{
			rho_star = state_j.rho_;
		}

		return rho_star;
	}
	//=================================================================================================//
	Real HLLCWithLimiterRiemannSolver::
		getEStar(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = dot(-e_ij, state_i.vel_);
		Real ur = dot(-e_ij, state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
		Real energy_star = 0.0;
		if (0.0 < s_l)
		{
			energy_star = state_i.E_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			energy_star = state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)));
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			energy_star = state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
		}
		if (s_r < 0.0)
		{
			energy_star = state_j.E_;
		}

		return energy_star;
	}
	//=================================================================================================//
}
