#ifndef RIEMANN_SOLVER_EXTRA_H
#define RIEMANN_SOLVER_EXTRA_H

#include "riemann_solver.h"

namespace SPH
{
	class AcousticRiemannSolverExtra : public AcousticRiemannSolver
	{
	public:
		template <class FluidI, class FluidJ>
		AcousticRiemannSolverExtra(FluidI& fluid_i, FluidJ& fluid_j)
			: AcousticRiemannSolver(fluid_i, fluid_j) {};
		Real DissipativePJump(const Real& u_jump);
	};

} // namespace SPH
#endif // RIEMANN_SOLVER_EXTRA_H
