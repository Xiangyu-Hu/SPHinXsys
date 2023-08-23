#include "riemann_solver_extra.h"

namespace SPH
{
	Vecd AcousticRiemannSolverExtra::DissipativePJumpExtra(const Vecd& u_jump, const Vecd& e_ij)
	{
		return rho0c0_geo_ave_ * u_jump * SMIN(Real(3) * SMAX(u_jump.dot(e_ij) * inv_c_ave_, Real(0)), Real(1));
	}
	Vecd DissipativeRiemannSolverExtra::DissipativePJumpExtra(const Vecd& u_jump, const Vecd& e_ij)
	{
		return rho0c0_geo_ave_ * u_jump;
	}
} // namespace SPH