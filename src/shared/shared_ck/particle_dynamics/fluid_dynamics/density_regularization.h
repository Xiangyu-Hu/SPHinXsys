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

template <typename... RelationTypes>
class DensityRegularization;

template <template <typename...> class RelationType, typename... Parameters>
class DensityRegularization<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class DynamicsIdentifier>
    explicit DensityRegularization(DynamicsIdentifier &identifier);
    virtual ~DensityRegularization() {};

    class InteractKernel : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       DensityRegularization<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);
        Real InitialDensity() { return rho0_; };

      protected:
        Real *rho_, *mass_, *rho_sum_, *Vol_;
        Real rho0_, inv_sigma0_;
    };

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_rho_sum_;
    Real rho0_, inv_sigma0_;
};

template <class FlowType, class ParticleScopeType, typename... Parameters>
class DensityRegularization<Inner<WithUpdate, FlowType, ParticleScopeType, Parameters...>>
    : public DensityRegularization<Base, Inner<Parameters...>>
{
    using RegularizationKernel = typename Regularization<FlowType>::ComputingKernel;
    using ParticleScopeTypeKernel = typename ParticleScopeTypeCK<ParticleScopeType>::ComputingKernel;

  public:
    explicit DensityRegularization(Inner<Parameters...> &inner_relation);
    virtual ~DensityRegularization() {};

    class InteractKernel
        : public DensityRegularization<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       DensityRegularization<Inner<WithUpdate, FlowType, ParticleScopeType, Parameters...>> &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real W0_;
    };

    class UpdateKernel
        : public DensityRegularization<Base, Inner<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        UpdateKernel(const ExecutionPolicy &ex_policy,
                     DensityRegularization<Inner<WithUpdate, FlowType, ParticleScopeType, Parameters...>> &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        RegularizationKernel regularization_;
        ParticleScopeTypeKernel particle_scope_;
    };

  protected:
    Regularization<FlowType> regularization_method_;
    ParticleScopeTypeCK<ParticleScopeType> within_scope_method_;
};

template <typename... Parameters>
class DensityRegularization<Contact<Parameters...>>
    : public DensityRegularization<Base, Contact<Parameters...>>
{
  public:
    explicit DensityRegularization(Contact<Parameters...> &contact_relation);
    virtual ~DensityRegularization() {};

    class InteractKernel
        : public DensityRegularization<Base, Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       DensityRegularization<Contact<Parameters...>> &encloser,
                       size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real contact_inv_rho0_k_;
        Real *contact_mass_k_;
    };

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<DiscreteVariable<Real> *> dv_contact_mass_;
};

using DensityRegularizationComplex = DensityRegularization<Inner<WithUpdate, Internal, AllParticles>, Contact<>>;
using DensityRegularizationComplexFreeSurface = DensityRegularization<Inner<WithUpdate, FreeSurface, AllParticles>, Contact<>>;
using DensityRegularizationComplexInternalPressureBoundary = DensityRegularization<Inner<WithUpdate, Internal, ExcludeBufferParticles>, Contact<>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // DENSITY_REGULARIZATION_H
