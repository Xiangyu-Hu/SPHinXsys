/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *                                                         *
 * -----------------------------------------------------------------------------*/

#ifndef GENERAL_DYNAMICS_WKGC_H
#define GENERAL_DYNAMICS_WKGC_H

#include "general_dynamics.h"

#include <limits>

namespace SPH
{
	typedef DataDelegateInner<BaseParticles> GeneralDataDelegateInner;
	typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase>
		GeneralDataDelegateContact;

	class GlobalCorrectionMatrix : public LocalDynamics, public GeneralDataDelegateInner, public GeneralDataDelegateContact
	{
	public:
		GlobalCorrectionMatrix(ComplexRelation& complex_relation, int beta, Real alpha)
			: LocalDynamics(complex_relation.getSPHBody()),
			GeneralDataDelegateInner(complex_relation.getInnerRelation()),
			GeneralDataDelegateContact(complex_relation.getContactRelation()),
			beta_(beta), alpha_(alpha), Vol_(particles_->Vol_)
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
		const int beta_;
		const Real alpha_;
		StdLargeVec<Real>& Vol_;
		StdLargeVec<Matd>B_;

		void interaction(size_t index_i, Real dt = 0.0);
	};
}
#endif // GENERAL_DYNAMICS_WKGC_H