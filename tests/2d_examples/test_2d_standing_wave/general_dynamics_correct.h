/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4	.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/

#ifndef GENERAL_DYNAMICS_CORRECT_H
#define GENERAL_DYNAMICS_CORRECT_H

#include "all_particle_dynamics.h"
#include "base_body.h"
#include "base_particles.h"
#include "external_force.h"

#include <limits>

namespace SPH
{
	typedef DataDelegateInner<BaseParticles> GeneralDataDelegateInner;
	typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase>
		GeneralDataDelegateContact;

	class GlobalCorrectionMatrix : public LocalDynamics, public GeneralDataDelegateInner, public GeneralDataDelegateContact
	{
	public:
		GlobalCorrectionMatrix(ComplexRelation& complex_relation)
			: LocalDynamics(complex_relation.getSPHBody()),
			GeneralDataDelegateInner(complex_relation.getInnerRelation()),
			GeneralDataDelegateContact(complex_relation.getContactRelation()),
			Vol_(particles_->Vol_)
		{
			particles_->registerVariable(B_, "WeightedCorrectionMatrix");

			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_mass_.push_back(&(contact_particles_[k]->mass_));
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}

		virtual ~GlobalCorrectionMatrix() {}

	protected:
		StdVec<StdLargeVec<Real>*> contact_Vol_;
		StdVec<StdLargeVec<Real>*> contact_mass_;
		StdLargeVec<Real>& Vol_;
		StdLargeVec<Matd>B_;

		void interaction(size_t index_i, Real dt = 0.0);
	};
}
#endif // GENERAL_DYNAMICS_H