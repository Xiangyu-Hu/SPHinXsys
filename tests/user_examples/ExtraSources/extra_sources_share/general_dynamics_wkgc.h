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
	typedef DataDelegateContact<BaseParticles, BaseParticles> GeneralDataDelegateFullContact;
class CorrectionMatrixInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    CorrectionMatrixInner(BaseInnerRelation &inner_relation, int beta, Real alpha)
        : LocalDynamics(inner_relation.getSPHBody()),
          GeneralDataDelegateInner(inner_relation),
          beta_(beta), alpha_(alpha),
          B_(*particles_->registerSharedVariable<Matd>("WeightedCorrectionMatrix")){};
    CorrectionMatrixInner(LocalDynamicsParameters<BaseInnerRelation, int, Real> parameter_set)
        : CorrectionMatrixInner(parameter_set.body_relation_,
                                std::get<0>(parameter_set.others_),
                                std::get<1>(parameter_set.others_)){};

    virtual ~CorrectionMatrixInner() {}

  protected:
    int beta_;
    Real alpha_;
    StdLargeVec<Matd> &B_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

class CorrectionMatrixContact : public LocalDynamics, public GeneralDataDelegateFullContact
{
  public:
    CorrectionMatrixContact(BaseContactRelation &contact_relation)
        : LocalDynamics(contact_relation.getSPHBody()),
          GeneralDataDelegateFullContact(contact_relation),
          B_(*particles_->getVariableByName<Matd>("WeightedCorrectionMatrix"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_mass_.push_back(&(contact_particles_[k]->mass_));
            contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
        }
    }

    virtual ~CorrectionMatrixContact() {}

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    StdLargeVec<Matd> &B_;

    void interaction(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // GENERAL_DYNAMICS_WKGC_H