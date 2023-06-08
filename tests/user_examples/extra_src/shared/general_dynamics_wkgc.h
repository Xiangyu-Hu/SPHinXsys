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
class CorrectionMatrixInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    CorrectionMatrixInner(BaseInnerRelation &inner_relation, int beta, Real alpha);
    virtual ~CorrectionMatrixInner(){};

  protected:
    int beta_;
    Real alpha_;
    StdLargeVec<Matd> &B_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

class CorrectionMatrixComplex : public CorrectionMatrixInner, public GeneralDataDelegateContact
{
  public:
    CorrectionMatrixComplex(ComplexRelation &complex_relation, int beta, Real alpha);
    virtual ~CorrectionMatrixComplex(){};

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;

    void interaction(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // GENERAL_DYNAMICS_WKGC_H