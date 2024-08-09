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
 * @file non_reflective_boundary.h
 * @brief Here, we define the non-reflective boundary condition at
 * fluid surface particles. //TODO: extend the method so that it works for
 * both eulerian and lagrangian formulations.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef NON_REFLECTIVE_BOUNDARY_H
#define NON_REFLECTIVE_BOUNDARY_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
class NonReflectiveBoundaryCorrection : public LocalDynamics, public DataDelegateInner
{
  public:
    NonReflectiveBoundaryCorrection(BaseInnerRelation &inner_relation);
    virtual ~NonReflectiveBoundaryCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
    Real rho_farfield_, sound_speed_;
    Vecd vel_farfield_;
    Real *rho_, *p_, *Vol_, *mass_;
    Vecd *vel_, *mom_, *pos_;
    Real *inner_weight_summation_, *rho_average_, *vel_normal_average_;
    Vecd *vel_tangential_average_, *vel_average_;
    int *indicator_, *smeared_surface_;
    Vecd *n_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // NON_REFLECTIVE_BOUNDARY_H
