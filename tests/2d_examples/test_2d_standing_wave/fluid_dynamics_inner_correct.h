/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/

#ifndef FLUID_DYNAMICS_INNER_CORRECT_H
#define FLUID_DYNAMICS_INNER_CORRECT_H

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "fluid_body.h"
#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"
#include "riemann_solver.h"
#include "fluid_dynamics_inner.h"

namespace SPH
{
	namespace fluid_dynamics
	{
		/**
		* @class BaseIntegration1stHalfCorrect
		* @brief Template class for pressure relaxation scheme with the Riemann solver
		* as template variable
		*/
		template <class RiemannSolverType>
		class BaseIntegration1stHalfCorrect : public BaseIntegration
		{
		public:
			explicit BaseIntegration1stHalfCorrect(BaseInnerRelation& inner_relation);
			virtual ~BaseIntegration1stHalfCorrect() {};
			RiemannSolverType riemann_solver_;
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			virtual Vecd computeNonConservativeAcceleration(size_t index_i);
			StdLargeVec<Matd> p_B_;
			StdLargeVec<Matd>& B_;
		};
		using Integration1stHalfCorrect = BaseIntegration1stHalfCorrect<NoRiemannSolver>;
		/** define the mostly used pressure relaxation scheme using Riemann solver */
		using Integration1stHalfRiemannCorrect = BaseIntegration1stHalfCorrect<AcousticRiemannSolver>;
		using Integration1stHalfDissipativeRiemannCorrect = BaseIntegration1stHalfCorrect<DissipativeRiemannSolver>;
	}
}
#endif // FLUID_DYNAMICS_INNER_H