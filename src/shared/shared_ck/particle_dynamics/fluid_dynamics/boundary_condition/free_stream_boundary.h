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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	free_stream_boundary.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef FREE_STREAM_BOUNDARY_H
#define FREE_STREAM_BOUNDARY_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename ConditionFunction>
class FreeStreamCondition : public LocalDynamics
{
  public:
    explicit FreeStreamCondition(SPHBody &sph_body, const ConditionFunction &free_stream_velocity);
    virtual ~FreeStreamCondition() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        ConditionFunction free_stream_velocity_;
        Real rho0_;
        Real *rho_sum_;
        Vecd *pos_, *vel_;
        int *indicator_;
        Real *physical_time_;
    };

  protected:
    ConditionFunction free_stream_velocity_;
    Real rho0_;
    DiscreteVariable<Real> *dv_rho_sum_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_;
    DiscreteVariable<int> *dv_indicator_;
    SingularVariable<Real> *sv_physical_time_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FREE_STREAM_BOUNDARY_H