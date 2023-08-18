/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	general_interaction.h
 * @brief 	This is the interaction dynamics applicable for all type bodies
 * @author	Yaru Ren and Xiangyu Hu
 */

#ifndef GENERAL_INTERACTION_H
#define GENERAL_INTERACTION_H

#include "general_dynamics.h"

namespace SPH
{
class CorrectedConfigurationInner : public LocalDynamics, public GeneralDataDelegateInner
{
  public:
    CorrectedConfigurationInner(BaseInnerRelation &inner_relation, int beta = 0, Real alpha = Real(0));
    virtual ~CorrectedConfigurationInner(){};

  protected:
    int beta_;
    Real alpha_;
    StdLargeVec<Matd> &B_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

class CorrectedConfigurationComplex : public CorrectedConfigurationInner, public GeneralDataDelegateContactOnly
{
  public:
    CorrectedConfigurationComplex(ComplexRelation &complex_relation, int beta = 0, Real alpha = Real(0));
    virtual ~CorrectedConfigurationComplex(){};

  protected:
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> contact_mass_;

    void interaction(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // GENERAL_INTERACTION_H