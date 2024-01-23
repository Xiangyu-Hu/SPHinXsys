#include "riemann_solver_extra.h"

namespace SPH
{
	Real AcousticRiemannSolverExtra::DissipativePJump(const Real& u_jump)
	{
		return rho0c0_geo_ave_ * u_jump * SMIN((Real)Dimensions * Real(20) * SMAX(u_jump * inv_c_ave_, Real(0)), Real(1));
	}
} // namespace SPH