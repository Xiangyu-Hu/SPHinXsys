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
/**
 * @file 	general_solid_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef GENERAL_SOLID_DYNAMICS_H
#define GENERAL_SOLID_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidParticles> SolidDataSimple;
		typedef DataDelegateInner<SolidParticles> SolidDataInner;

		/**
		 * @class CorrectConfiguration
		 * @brief obtain the corrected initial configuration in strong form
		 */
		class CorrectConfiguration : public LocalDynamics, public SolidDataInner
		{
		public:
			explicit CorrectConfiguration(BaseInnerRelation &inner_relation);
			virtual ~CorrectConfiguration(){};

			inline void interaction(size_t index_i, Real dt = 0.0)
			{
				Matd local_configuration = Eps * Matd::Identity(); // a small number added to diagonal to avoid divide zero
				const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
					Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
					local_configuration -= r_ji * gradW_ijV_j.transpose();
				}
				B_[index_i] = local_configuration.inverse();
			};

		protected:
			StdLargeVec<Matd> &B_;
		};
	}
}
#endif // GENERAL_SOLID_DYNAMICS_H