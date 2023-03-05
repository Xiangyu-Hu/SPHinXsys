#include "riemann_solver.h"
#include "base_material.h"
#include "compressible_fluid.h"

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
		return rho0c0_geo_ave_ * u_jump * SMIN(3.0 * SMAX(u_jump * inv_c_ave_, 0.0), 1.0);
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
	FluidState HLLCRiemannSolverInWeaklyCompressibleFluid::
		getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		Real ul = -e_ij.dot(state_i.vel_);
		Real ur = -e_ij.dot(state_j.vel_);
		Real s_l = ul - fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l) + TinyReal);
		Real p_star = 0.0;
		Vecd v_star = Vecd::Zero();
		Real rho_star = 0.0;
		if (0.0 < s_l)
		{
			p_star = state_i.p_;
			v_star = state_i.vel_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
			v_star = state_i.vel_ - e_ij * (s_star - ul);
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
			v_star = state_j.vel_ - e_ij * (s_star - ur);
		}
		if (s_r < 0.0)
		{
			p_star = state_j.p_;
			v_star = state_j.vel_;
		}

		FluidState interface_state(rho_star, v_star, p_star);
		interface_state.vel_ = v_star;
		interface_state.p_ = p_star;

		return interface_state;
	}
	//=================================================================================================//
	FluidStarState HLLCRiemannSolverWithLimiterInWeaklyCompressibleFluid::
		getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij)
	{
		Real ul = -e_ij.dot(state_i.vel_);
		Real ur = -e_ij.dot(state_j.vel_);
		Real s_l = ul - fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l) + TinyReal);
		Real p_star = 0.0;
		Vecd v_star = Vecd::Zero();
		if (0.0 < s_l)
		{
			p_star = state_i.p_;
			v_star = state_i.vel_;
		}
		if (s_l <= 0.0 && 0.0 <= s_star)
		{
			Real rho_ave = 2 * state_i.rho_ * state_j.rho_ / (state_i.rho_ + state_j.rho_);
			Real rho_cl = state_i.rho_ * fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
			Real rho_cr = state_j.rho_ * fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
			Real rho_clr = (rho_cl * state_i.rho_ + rho_cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);
			p_star = 0.5 * (state_i.p_ + state_j.p_) + 0.5 * (SMIN(3.0 * SMAX(rho_ave * (ul - ur), 0.0), rho_clr) * (ul - ur) + s_star * (rho_cr - rho_cl));
			v_star = state_i.vel_ - e_ij * (s_star - ul);
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			Real rho_ave = 2 * state_i.rho_ * state_j.rho_ / (state_i.rho_ + state_j.rho_);
			Real rho_cl = state_i.rho_ * fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
			Real rho_cr = state_j.rho_ * fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
			Real rho_clr = (rho_cl * state_i.rho_ + rho_cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);
			p_star = 0.5 * (state_i.p_ + state_j.p_) + 0.5 * (SMIN(3.0 * SMAX(rho_ave * (ul - ur), 0.0), rho_clr) * (ul - ur) + s_star * (rho_cr - rho_cl));
			v_star = state_j.vel_ - e_ij * (s_star - ur);
		}
		if (s_r < 0.0)
		{
			p_star = state_j.p_;
			v_star = state_j.vel_;
		}

		FluidStarState interface_state(v_star, p_star);

		return interface_state;
	}
	//=================================================================================================//
	CompressibleFluidStarState HLLCRiemannSolver::
		getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = -e_ij.dot(state_i.vel_);
		Real ur = -e_ij.dot(state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
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
			energy_star = state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)));
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
			v_star = state_j.vel_ - e_ij * (s_star - ur);
			rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
			energy_star = state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
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
	CompressibleFluidStarState HLLCWithLimiterRiemannSolver::
		getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
	{
		Real ul = -e_ij.dot(state_i.vel_);
		Real ur = -e_ij.dot(state_j.vel_);
		Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
		Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
		Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
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
			Real rho_ave = 2 * state_i.rho_ * state_j.rho_ / (state_i.rho_ + state_j.rho_);
			Real rho_cl = state_i.rho_ * compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
			Real rho_cr = state_j.rho_ * compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
			Real rho_clr = (rho_cl * state_i.rho_ + rho_cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);
			p_star = 0.5 * (state_i.p_ + state_j.p_) + 0.5 * (SMIN(3.0 * SMAX(rho_ave * (ul - ur), 0.0), rho_clr) * (ul - ur) + s_star * (rho_cr - rho_cl));
			v_star = state_i.vel_ - e_ij * (s_star - ul);
			rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
			energy_star = state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)));
		}
		if (s_star <= 0.0 && 0.0 <= s_r)
		{
			Real rho_ave = 2 * state_i.rho_ * state_j.rho_ / (state_i.rho_ + state_j.rho_);
			Real rho_cl = state_i.rho_ * compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
			Real rho_cr = state_j.rho_ * compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
			Real rho_clr = (rho_cl * state_i.rho_ + rho_cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);
			p_star = 0.5 * (state_i.p_ + state_j.p_) + 0.5 * (SMIN(3.0 * SMAX(rho_ave * (ul - ur), 0.0), rho_clr) * (ul - ur) + s_star * (rho_cr - rho_cl));
			v_star = state_j.vel_ - e_ij * (s_star - ur);
			rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
			energy_star = state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
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
}
