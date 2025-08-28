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
 * @file force_prior_ck.h
 * @brief Here, we define the base methods for force prior
 * (all forces on particle except that due to pressure or stress)
 * used in SPH method.
 * @author Xiangyu Hu
 */
#ifndef FORCE_PRIOR_CK_H
#define FORCE_PRIOR_CK_H

#include "base_general_dynamics.h"

namespace SPH
{

class ForcePriorCK
{
  public:
    ForcePriorCK(BaseParticles *particles, const std::string &force_name)
        : dv_force_prior_(particles->registerStateVariableOnly<Vecd>("ForcePrior")),
          dv_current_force_(particles->registerStateVariableOnly<Vecd>(force_name)),
          dv_previous_force_(particles->registerStateVariableOnly<Vecd>("Previous" + force_name))
    {
        particles->addEvolvingVariable<Vecd>("Previous" + force_name);
    };
    virtual ~ForcePriorCK() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0)
        {
            force_prior_[index_i] += current_force_[index_i] - previous_force_[index_i];
            previous_force_[index_i] = current_force_[index_i];
        };

      protected:
        Vecd *force_prior_, *current_force_, *previous_force_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_force_prior_, *dv_current_force_, *dv_previous_force_;
};

template <class GravityType>
class GravityForceCK : public LocalDynamics, public ForcePriorCK
{
  public:
    GravityForceCK(SPHBody &sph_body, const GravityType &gravity);
    virtual ~GravityForceCK() {};

    class UpdateKernel : public ForcePriorCK::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        GravityType gravity_;
        Real *physical_time_;
        Vecd *pos_;
        Real *mass_;
    };

  protected:
    const GravityType gravity_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Real> *dv_mass_;
};
} // namespace SPH
#endif // FORCE_PRIOR_CK_H
