/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                        *
 * ------------------------------------------------------------------------*/

#ifndef FLUID_DYNAMICS_INNER_WKGC_H
#define FLUID_DYNAMICS_INNER_WKGC_H

#include "fluid_dynamics_inner.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		/**
		* @class BaseIntegration1stHalfCorrect
		*/
		template <class RiemannSolverType>
		class BaseIntegration1stHalfCorrect : public BaseIntegration1stHalf<RiemannSolverType>
		{
		public:
			explicit BaseIntegration1stHalfCorrect(BaseInnerRelation& inner_relation);
			virtual ~BaseIntegration1stHalfCorrect() {};

			using BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf;
			void initialization(size_t index_i, Real dt);
			void interaction(size_t index_i, Real dt);

		protected:
			StdLargeVec<Matd> p_B_;
			StdLargeVec<Matd>& B_;
		};
		using Integration1stHalfCorrect = BaseIntegration1stHalfCorrect<NoRiemannSolver>;
		/** define the mostly used pressure relaxation scheme using Riemann solver */
		using Integration1stHalfRiemannCorrect = BaseIntegration1stHalfCorrect<AcousticRiemannSolver>;
	}
}
#endif // FLUID_DYNAMICS_INNER_WKGC_H