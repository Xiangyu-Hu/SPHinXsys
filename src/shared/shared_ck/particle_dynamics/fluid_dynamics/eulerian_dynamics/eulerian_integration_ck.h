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
 * @file 	eulerian_integration_ck.h
 * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_INTEGRATION_CK_H
#define EULERIAN_INTEGRATION_CK_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class BaseInteractionType>
class EulerianIntegrationCK : public BaseInteractionType
{
  public:
    template <class BaseRelationType>
    explicit EulerianIntegrationCK(BaseRelationType &base_relation);
    virtual ~EulerianIntegrationCK() {};

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_p_, *dv_dmass_dt_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_mom_, *dv_dmom_dt_;
};

template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
class EulerianIntegrationCK<Inner<RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public EulerianIntegrationCK<Inner<Parameters...>>
{
    using EosKernel = typename WeaklyCompressibleFluid::EosKernel;
    using BaseInteraction = EulerianIntegrationCK<Inner<Parameters...>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    template <class BaseRelationType>
    explicit EulerianIntegrationCK(BaseRelationType &base_relation);
    virtual ~EulerianIntegrationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *rho_, *p_, *Vol_, *dmass_dt_;
        Vecd *vel_, *mom_, *dmom_dt_;
    };

  protected:
    KernelCorrectionType kernel_correction_method_;
    WeaklyCompressibleFluid &fluid_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, class... Parameters>
class EulerianIntegrationCK<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public EulerianIntegrationCK<Contact<Parameters...>>, public Interaction<Wall>
{
    using BaseInteraction = EulerianIntegrationCK<Contact<Parameters...>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    template <class BaseRelationType>
    explicit EulerianIntegrationCK(BaseRelationType &base_relation);
    virtual ~EulerianIntegrationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *rho_, *p_, *Vol_, *dmass_dt_;
        Vecd *vel_, *mom_, *dmom_dt_;
        Vecd *wall_n_, *wall_vel_ave_;
    };

  protected:
    KernelCorrectionType kernel_correction_method_;
    RiemannSolverType riemann_solver_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
    : public EulerianIntegrationCK<RelationType<InteractionParameters...>>
{
    using EosKernel = typename WeaklyCompressibleFluid::EosKernel;
    using BaseDynamicsType = EulerianIntegrationCK<RelationType<InteractionParameters...>>;

  public:
    template <class BaseRelationType, typename... Args>
    explicit EulerianIntegrationCK(BaseRelationType &base_relation, Args &&...args);
    virtual ~EulerianIntegrationCK() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        Real *dmass_dt_;
        Vecd *dmom_dt_;
    };

    class UpdateKernel : public BaseInteraction::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0);

      protected:
        EosKernel eos_;
        Real *rho_, *p_, *Vol_, *dmass_dt_;
        Vecd *vel_, *mom_, *dmom_dt_, *force_prior_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_force_prior_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class EulerianIntegrationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>
    : public EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
{
    using BaseDynamicsType = EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>;

  public:
    template <typename... Args>
    explicit EulerianIntegrationCK(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~EulerianIntegrationCK() {};
};

template <template <typename...> class RelationType, class... InteractionParameters>
class EulerianIntegrationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>
    : public EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
{
    using BaseDynamicsType = EulerianIntegrationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>;

  public:
    template <typename... Args>
    explicit EulerianIntegrationCK(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~EulerianIntegrationCK() {};
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_INTEGRATION_CK_H
