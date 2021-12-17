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
* @file 	complex_solid.h
* @brief 	These are classes for define complex solid materials.
* @author	Xiangyu Hu and Chi Zhang
*/
#pragma once

#include "elastic_solid.h"

namespace SPH
{

	class ActiveMuscleParticles;

	/**
	* @class ActiveMuscle
	* @brief Here, the active reponse is considered.
	*/
	template <class MuscleType>
	class ActiveMuscle : public MuscleType
	{
	protected:
		ActiveMuscleParticles *active_muscle_particles_;

	public:
		template <typename... ConstructorArgs>
		explicit ActiveMuscle(ConstructorArgs &&...args);
		virtual ~ActiveMuscle(){};

		void assignActiveMuscleParticles(ActiveMuscleParticles *active_muscle_particles);
		/** compute the stress through Constitutive relation. */
		virtual Matd ConstitutiveRelation(Matd &deformation, size_t index_i) override;
		virtual ActiveMuscle<MuscleType> *ThisObjectPtr() override { return this; };
	};
}
