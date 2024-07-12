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
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	smeared_surface_indication.h
 * @brief 	Here, we define methods for smear material surface,
 * i.e. the surface is smeared to the thickness of a cutoff radius.
 * @author	Xiangyu Hu
 */
#ifndef SMEARED_SURFACE_INDICATION_H
#define SMEARED_SURFACE_INDICATION_H

#include "base_general_dynamics.h"

namespace SPH
{
class SmearedSurfaceIndication : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit SmearedSurfaceIndication(BaseInnerRelation &inner_relation);
    virtual ~SmearedSurfaceIndication(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    int *indicator_;
    int *smeared_surface_;
};
} // namespace SPH
#endif // SMEARED_SURFACE_INDICATION_H