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
 *  HU1527/12-1 and HU1527/12-4												*
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
 * @file 	complex_solid.h
 * @brief 	These are classes for define complex solid materials. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#pragma once

#include "elastic_solid.h"

namespace SPH
{
	/**
	 * @class ActiveMuscle
	 * @brief Here, the active response is considered.
	 */
	template <class MuscleType>
	class ActiveMuscle : public MuscleType
	{
	protected:
		StdLargeVec<Real> active_contraction_stress_; /**<  active contraction stress */

	public:
		template <typename... ConstructorArgs>
		explicit ActiveMuscle(ConstructorArgs &&...args);
		virtual ~ActiveMuscle(){};

		/** initialize the local properties, fiber and sheet direction. */
		virtual void initializeLocalParameters(BaseParticles *base_particles) override;
		/** compute the stress through Constitutive relation. */
		virtual Matd StressPK2(Matd &deformation, size_t index_i) override;
		virtual ActiveMuscle<MuscleType> *ThisObjectPtr() override { return this; };
	};
}
