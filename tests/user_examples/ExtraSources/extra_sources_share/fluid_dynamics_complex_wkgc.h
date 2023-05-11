/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*																		*                                                                        *
 * ------------------------------------------------------------------------*/

#ifndef FLUID_DYNAMICS_COMPLEX_WKGC_H
#define FLUID_DYNAMICS_COMPLEX_WKGC_H

#include "fluid_dynamics_inner_wkgc.h"
#include "fluid_dynamics_inner_wkgc.hpp"
#include "fluid_dynamics_complex.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		/**
        * @class BaseIntegration1stHalfCorrectWithWall
        * @brief  template class pressure relaxation scheme together with wall boundary
        */
		template <class BaseIntegration1stHalfCorrectType>
		class BaseIntegration1stHalfCorrectWithWall : public InteractionWithWall<BaseIntegration1stHalfCorrectType>
		{
		public:
			template <typename... Args>
			BaseIntegration1stHalfCorrectWithWall(Args &&...args)
				: InteractionWithWall<BaseIntegration1stHalfCorrectType>(std::forward<Args>(args)...) {};
			virtual ~BaseIntegration1stHalfCorrectWithWall() {};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		using Integration1stHalfCorrectWithWall = BaseIntegration1stHalfCorrectWithWall<Integration1stHalfCorrect>;
		using Integration1stHalfRiemannCorrectWithWall = BaseIntegration1stHalfCorrectWithWall<Integration1stHalfRiemannCorrect>;
	}
}
#endif // FLUID_DYNAMICS_COMPLEX_WKGC_H