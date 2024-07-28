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
 * @file force_prior.h
 * @brief Here, we define the base methods for force prior
 * (all forces on particle except that due to pressure or stress)
 * used in SPH method.
 * @author Xiangyu Hu
 */
#ifndef FORCE_PRIOR_H
#define FORCE_PRIOR_H

#include "base_general_dynamics.h"

namespace SPH
{

class ForcePrior
{
  protected:
    Vecd *force_prior_, *current_force_, *previous_force_;

    template <>
    class ComputingKernel<ForcePrior>
    {
        Vecd *force_prior_, *current_force_, *previous_force_;

      public:
        ComputingKernel(BaseParticles *base_particles, const std::string &force_name);
        void update(size_t index_i, Real dt = 0.0);
    };

  public:
    ForcePrior(BaseParticles *base_particles, const std::string &force_name);
    template <class ExecutionPolicy>
    ForcePrior(const ExecutionPolicy &device_policy, BaseParticles *base_particles, const std::string &force_name);
    virtual ~ForcePrior(){};
    void update(size_t index_i, Real dt = 0.0);

    ComputingKernel<ForcePrior> *getComputingKernel() { return &computing_kernel_; }

  protected:
    ComputingKernel<ForcePrior> computing_kernel_;
};

template <class GravityType>
class GravityForce : public LocalDynamics, public DataDelegateSimple, public ForcePrior
{
  protected:
    GravityType gravity_;
    Vecd *pos_;
    Real *mass_;
    Real *physical_time_;

  public:
    explicit GravityForce(SPHBody &sph_body, const GravityType &gravity);
    virtual ~GravityForce(){};
    void update(size_t index_i, Real dt = 0.0);
};

} // namespace SPH
#endif // FORCE_PRIOR_H
