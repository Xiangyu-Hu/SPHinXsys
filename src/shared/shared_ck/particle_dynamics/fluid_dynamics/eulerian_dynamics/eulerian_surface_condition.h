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
 * @file eulerian_surface_condition.h
 * @brief TBD.
 * @author Xiangyu Hu
 */
#ifndef EULERIAN_SURFACE_CONDITION_H
#define EULERIAN_SURFACE_CONDITION_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class KernelCorrectionType, typename ConditionType>
class EulerianSurfaceCondition : public BaseLocalDynamics<BodyPartByParticle>,
                                 public EulerianIntegrationCK
{
    using EosKernel = typename WeaklyCompressibleFluid::EosKernel;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    template <typename... Args>
    EulerianSurfaceCondition(BodyPartByParticle &box_part, Args &&...args);

    class UpdateKernel : public EulerianIntegrationCK::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_kernel_;
        EosKernel eos_;
        ConditionType condition_;
        Real *physical_time_;
        Vecd *n_, *vel_, *dmom_dt_;
        Real *surface_area_, *rho_, *p_, *dmass_dt_;
    };

  protected:
    KernelCorrectionType kernel_correction_method_;
    WeaklyCompressibleFluid &fluid_;
    ConditionType condition_;
    SingularVariable<Real> *sv_physical_time_;
    DiscreteVariable<Vecd> *dv_n_;
    DiscreteVariable<Real> *dv_surface_area_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_SURFACE_CONDITION_H
