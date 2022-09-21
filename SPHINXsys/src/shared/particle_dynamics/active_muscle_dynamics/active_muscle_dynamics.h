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
 * @file 	active_muscle_dynamics.h
 * @brief 	In is file, we declare muscle dynamics which is driven by an external injection of energy. 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.3.1
 *			Here, we need identify the physical differences between electrophysiology and active muscle.
 *			The former is on the diffusion and electro-chemical reaction happens in tissue.
 *			The latter is for muscle dynamics which is driven by an external injection of energy.
 *			As the naming of class, function and variables in this code need be based on physics.
 *			We will identify physical differences by properly choosing names.
 *			Xiangyu Hu
 */

#ifndef ACTIVE_MUSCLE_DYNAMICS_H
#define ACTIVE_MUSCLE_DYNAMICS_H



#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "base_kernel.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
    namespace active_muscle_dynamics
	{
        typedef DataDelegateSimple<ElasticSolidParticles> ElasticSolidDataSimple;
 
		/**
		 * @class MuscleActivation
		 * @brief  impose cases specific muscle activation
		 * This is a abstract class to be override for case specific activation
		 */
		class MuscleActivation :
			public LocalDynamics, public ElasticSolidDataSimple
		{
		public:
			explicit MuscleActivation(SPHBody &sph_body);
			virtual ~MuscleActivation() {};
		protected:
			StdLargeVec<Vecd>& pos0_;
			StdLargeVec<Real>& active_contraction_stress_;
		};
    }
}
#endif //ACTIVE_MUSCLE_DYNAMICS_H