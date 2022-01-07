/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	inelastic_solid_dynamics.h
* @brief 	Here, we define the algorithm classes for inelastic_solid dynamics.
* @details 	We consider here a weakly compressible solids.
* @author	Xiaojing Tang, Chi Zhang and Xiangyu Hu
*/
#pragma once

#include "solid_dynamics.h"

namespace SPH
{
	namespace solid_dynamics
	{
		/**
		 * @class PlasticStressRelaxationFirstHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class PlasticStressRelaxationFirstHalf
			: public StressRelaxationFirstHalf
		{
		public:
			PlasticStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~PlasticStressRelaxationFirstHalf(){};

		protected:
			PlasticSolid *plastic_solid_;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
		};
	}
}
