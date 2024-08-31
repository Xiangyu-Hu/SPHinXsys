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
 * @file density_summation_ck.h
 * @brief Here, we define the algorithm classes for computing
 * the density of a continuum by kernel function summation.
 * @details We are using templates and their explicit or partial specializations
 * to identify variations of the interaction types..
 * @author Xiangyu Hu
 */

#ifndef DENSITY_SUMMATION_CK_H
#define DENSITY_SUMMATION_CK_H

#include "base_fluid_dynamics.h"
#include "local_interaction_dynamics_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... T>
class Regularization;

template <>
class Regularization<>
{
  public:
    template <class ComputingKernelType>
    Regularization(ComputingKernelType &computing_kernel){};
    Real operator()(Real &rho_sum, Real rho0, Real &rho) { return rho_sum; };
};

template <>
class Regularization<FreeSurface>
{
  public:
    template <class ComputingKernelType>
    Regularization(ComputingKernelType &computing_kernel){};
    Real operator()(Real &rho_sum, Real rho0, Real &rho) { return SMAX(rho_sum, rho0); };
};

template <typename... InteractionTypes>
class DensitySummationCK;

template <class LocalDynamicsType>
class DensitySummationCK<LocalDynamicsType> : public LocalDynamicsType
{
  public:
    template <class DynamicsIdentifier>
    explicit DensitySummationCK(DynamicsIdentifier &identifier);
    virtual ~DensitySummationCK(){};

    class ComputingKernel : public LocalDynamicsType::ComputingKernel
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        DensitySummationCK<LocalDynamicsType> &encloser);
        Real *rho_, *mass_, *rho_sum_, *Vol_;
        Real rho0_, inv_sigma0_;
    };

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_rho_sum_, *dv_Vol_;
    Real rho0_, inv_sigma0_;
};

template <typename... RegularizationType, typename... Parameters>
class DensitySummationCK<Inner<Regularization<RegularizationType...>, Parameters...>>
    : public DensitySummationCK<LocalInteractionDynamics<Inner<Parameters...>>>
{
  public:
    explicit DensitySummationCK(InnerRelation &inner_relation);
    virtual ~DensitySummationCK(){};

    class ComputingKernel
        : public DensitySummationCK<LocalInteractionDynamics<Inner<Parameters...>>>::
              ComputingKernel
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        DensitySummationCK<Inner<Parameters...>> &encloser);
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real W0_;
        Regularization<RegularizationType...> regularization_;
    };
};
using DensitySummationCKInner = DensitySummationCK<Inner<Regularization<>>>;

template <typename... Parameters>
class DensitySummationCK<Contact<Parameters...>>
    : public DensitySummationCK<LocalInteractionDynamics<Contact<Parameters...>>>
{
  public:
    explicit DensitySummationCK(ContactRelation &contact_relation);
    virtual ~DensitySummationCK(){};

    class ComputingKernel
        : public DensitySummationCK<LocalInteractionDynamics<Contact<Parameters...>>>::
              ComputingKernel
    {
      public:
        template <class ExecutionPolicy>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        DensitySummationCK<Contact<Parameters...>> &encloser,
                        size_t contact_index);
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        StdVec<Real> contact_inv_rho0_k_;
        StdVec<Real *> contact_mass_k_;
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<DiscreteVariable<Real> *> dv_contact_mass_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_CK_H
