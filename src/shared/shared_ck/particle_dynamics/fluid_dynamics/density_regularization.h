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
class CompressionSummation;

template <template <typename...> class RelationType, typename... Parameters>
class CompressionSummation<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class DynamicsIdentifier>
    explicit CompressionSummation(DynamicsIdentifier &identifier);
    virtual ~CompressionSummation() {};

  protected:
    DiscreteVariable<Real> *dv_Vol_ref_, *dv_compression_sum_;
};

template <typename... Parameters>
class CompressionSummation<Inner<Parameters...>>
    : public CompressionSummation<Base, Inner<Parameters...>>
{
    using BaseInteraction = CompressionSummation<Base, Inner<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit CompressionSummation(DynamicsIdentifier &identifier);
    virtual ~CompressionSummation() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        InteractKernel(const ExecutionPolicy &ex_policy, Encloser &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Vecd zero_;
        DataView<Real> Vol_ref_, compression_sum_;
    };
};

template <typename... Parameters>
class CompressionSummation<Contact<Parameters...>>
    : public CompressionSummation<Base, Contact<Parameters...>>
{
    using BaseInteraction = CompressionSummation<Base, Contact<Parameters...>>;

  public:
    template <class DynamicsIdentifier>
    explicit CompressionSummation(DynamicsIdentifier &identifier);
    virtual ~CompressionSummation() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class Encloser>
        InteractKernel(const ExecutionPolicy &ex_policy, Encloser &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        DataView<Real> compression_sum_, contact_Vol_ref_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_ref_;
};
//------------------------------------------------------------------
// forward declarations of Regularization<FlowType>
//------------------------------------------------------------------
template <typename...>
class Regularization;

template <class DynamicsIdentifier, class FluidType, class FlowType, typename... ParticleScopes>
class DensityRegularization : public BaseLocalDynamics<DynamicsIdentifier>
{
    using EosKernel = typename FluidType::EosKernel;
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
        void update(size_t index_i, Real dt = 0.0);

      protected:
        EosKernel eos_;
        DataView<Real> rho_, compression_sum_;
        RegularizationKernel regularization_;
        ParticleScopeTypeKernel particle_scope_;
    };

  protected:
    FluidType &fluid_;
    DiscreteVariable<Real> *dv_rho_, *dv_compression_sum_;
    Regularization<FlowType> regularization_method_;
    ParticleScopeTypeCK<ParticleScopes...> within_scope_method_;
};

template <>
class Regularization<Internal>
{
  public:
    Regularization(BaseParticles *base_particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser){};
        Real operator()(UnsignedInt index_i, const Real &compression_sum, const Real &rho0)
        {
            return compression_sum * rho0;
        };
    };
};

template <>
class Regularization<FreeSurface>
{
  public:
    Regularization(BaseParticles *base_particles) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser){};

        Real operator()(UnsignedInt index_i, const Real &compression_sum, const Real &rho0)
        {
            return SMAX(compression_sum, Real(1)) * rho0;
        };

      protected:
        Real rho0_;
    };
};

template <>
class Regularization<FreeStream>
{
    DiscreteVariable<int> *dv_indicator_;

  public:
    Regularization(BaseParticles *base_particles)
        : dv_indicator_(base_particles->getVariableByName<int>("Indicator")) {};

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser)
            : indicator_(encloser.dv_indicator_->DelegatedDataView(ex_policy)){};

        Real operator()(UnsignedInt index_i, const Real &compression_sum, const Real &rho0)
        {
            return indicator_[index_i] != 0 ? rho0 : compression_sum * rho0;
        };

      protected:
        DataView<int> indicator_;
    };
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_REGULARIZATION_H
