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
 * @file density_regularization.h
 * @brief Here, we define the algorithm classes for computing
 * the density of a continuum by kernel function summation.
 * @details We are using templates and their explicit or partial specializations
 * to identify variations of the interaction types..
 * @author Xiangyu Hu
 */

#ifndef DENSITY_REGULARIZATION_H
#define DENSITY_REGULARIZATION_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "particle_functors_ck.h" // or wherever ParticleScopeTypeCK is defined

namespace SPH
{
namespace fluid_dynamics
{
template <typename...>
class DensitySummationCK;

template <template <typename...> class RelationType, typename... Parameters>
class DensitySummationCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class DynamicsIdentifier>
    explicit DensitySummationCK(DynamicsIdentifier &identifier);
    virtual ~DensitySummationCK() {};

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class Encloser, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, Encloser &encloser, Args &&...args);
        Real InitialDensity() { return rho0_; };

      protected:
        Real *rho_, *mass_, *rho_sum_, *Vol_;
        Real rho0_;
    };

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_rho_sum_;
    Real rho0_;
};

template <typename... Parameters>
class DensitySummationCK<Inner<Parameters...>>
    : public DensitySummationCK<Base, Inner<Parameters...>>
{
  public:
    explicit DensitySummationCK(Inner<Parameters...> &inner_relation);
    virtual ~DensitySummationCK() {};

    class InteractKernel
        : public DensitySummationCK<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        InteractKernel(const ExecutionPolicy &ex_policy, Encloser &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Vecd zero_;
    };
};

template <typename... Parameters>
class DensitySummationCK<Contact<Parameters...>>
    : public DensitySummationCK<Base, Contact<Parameters...>>
{
  public:
    explicit DensitySummationCK(Contact<Parameters...> &contact_relation);
    virtual ~DensitySummationCK() {};

    class InteractKernel
        : public DensitySummationCK<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        InteractKernel(const ExecutionPolicy &ex_policy, Encloser &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real contact_inv_rho0_k_;
        Real *contact_mass_k_;
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<DiscreteVariable<Real> *> dv_contact_mass_;
};
//------------------------------------------------------------------
// forward declarations of Regularization<FlowType>
//------------------------------------------------------------------
template <typename...>
class Regularization;

template <>
class Regularization<Internal>
{
  public:
    Regularization(BaseParticles *particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        Regularization<Internal> &encloser,
                        ComputingKernelType &computing_kernel){};

        Real operator()(UnsignedInt index_i, Real &rho_sum) { return rho_sum; };
    };
};

template <>
class Regularization<FreeSurface>
{
  public:
    Regularization(BaseParticles *particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy,
                        Regularization<FreeSurface> &encloser,
                        ComputingKernelType &computing_kernel)
            : rho0_(computing_kernel.InitialDensity()){};

        Real operator()(UnsignedInt index_i, Real &rho_sum) { return SMAX(rho_sum, rho0_); };

      protected:
        Real rho0_;
    };
};

template <class DynamicsIdentifier, class FlowType, typename... ParticleScopes>
class DensityRegularization : public BaseLocalDynamics<DynamicsIdentifier>
{
    using RegularizationKernel = typename Regularization<FlowType>::ComputingKernel;
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopes...>::ComputingKernel;

  public:
    explicit DensityRegularization(DynamicsIdentifier &identifier);
    virtual ~DensityRegularization() {};

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        UpdateKernel(const ExecutionPolicy &ex_policy, Encloser &encloser);
        Real InitialDensity() { return rho0_; };
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real rho0_;
        Real *rho_, *rho_sum_;
        RegularizationKernel regularization_;
        ParticleScopeTypeKernel particle_scope_;
    };

  protected:
    Real rho0_;
    DiscreteVariable<Real> *dv_rho_, *dv_rho_sum_;
    Regularization<FlowType> regularization_method_;
    ParticleScopeTypeCK<ParticleScopes...> within_scope_method_;
};
} // namespace fluid_dynamics
} // namespace SPH

#endif // DENSITY_REGULARIZATION_H
