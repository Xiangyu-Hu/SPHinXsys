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
		Vecd DissipativePJumpExtra(const Vecd& u_jump, const Vecd& e_ij);
	};

	class DissipativeRiemannSolverExtra : public DissipativeRiemannSolver
	{
	public:
		template <class FluidI, class FluidJ>
		DissipativeRiemannSolverExtra(FluidI& fluid_i, FluidJ& fluid_j)
			: DissipativeRiemannSolver(fluid_i, fluid_j) {};
		Vecd DissipativePJumpExtra(const Vecd& u_jump, const Vecd& e_ij);
	};
} // namespace SPH
#endif // RIEMANN_SOLVER_EXTRA_H
