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
 * @file 	structure_dynamics.h
 * @brief 	TBD.
 * @author	Xiangyu Hu
 */

#ifndef STRUCTURE_DYNAMICS_H
#define STRUCTURE_DYNAMICS_H

#include "elastic_solid.h"
#include "force_prior_ck.h"
#include "interaction_ck.hpp"

namespace SPH
{
namespace solid_dynamics
{

class AcousticTimeStepCK : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit AcousticTimeStepCK(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStepCK() {};

    class FinishDynamics
    {
        Real h_min_, acousticCFL_;

      public:
        using OutputType = Real;
        FinishDynamics(AcousticTimeStepCK &encloser);
        Real Result(Real reduced_value);
    };

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Real reduce(size_t index_i, Real dt = 0.0);

      protected:
        Real h_min_, c0_;
        Real *mass_;
        Vecd *vel_, *force_, *force_prior_;
    };

  protected:
    Real acousticCFL_;
    Real h_min_, c0_;
    DiscreteVariable<Real> *dv_mass_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_force_, *dv_force_prior_;
};

class StructureDynamicsVariables
{
  public:
    explicit StructureDynamicsVariables(BaseParticles *particles);
    virtual ~StructureDynamicsVariables() {};

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_vel_, *dv_force_;
    DiscreteVariable<Matd> *dv_B_, *dv_F_, *dv_dF_dt_, *dv_inverse_F_, *dv_stress_on_particle_,
        *dv_scaling_matrix_;
};

template <typename...>
class StructureIntegration1stHalf;

template <class MaterialType, typename KernelCorrectionType, typename... Parameters>
class StructureIntegration1stHalf<Inner<OneLevel, MaterialType, KernelCorrectionType, Parameters...>>
    : public Interaction<Inner<Parameters...>>, public StructureDynamicsVariables
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;
    using Adaptation = typename Inner<Parameters...>::SourceType::Adaptation;
    using SmoothingLengthRatioType = typename Adaptation::SmoothingLengthRatioType;
    using ConstituteKernel = typename MaterialType::ConstituteKernel;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit StructureIntegration1stHalf(Inner<Parameters...> &inner_relation, Real numerical_damping_factor = 0.125);
    virtual ~StructureIntegration1stHalf() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        ConstituteKernel constitute_;
        CorrectionKernel correction_;
        Real rho0_, G_, h_ref_, numerical_damping_factor_;
        SmoothingLengthRatioType h_ratio_;
        Real *rho_;
        Vecd *pos_, *vel_;
        Matd *F_, *inverse_F_, *dF_dt_, *scaling_matrix_, *stress_on_particle_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real G_;
        Real *Vol0_;
        Vecd *pos_, *force_;
        Matd *scaling_matrix_, *inverse_F_, *stress_on_particle_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *mass_;
        Vecd *vel_, *force_, *force_prior_;
    };

  protected:
    MaterialType &material_;
    Adaptation &adaptation_;
    KernelCorrectionType kernel_correction_;
    Real h_ref_, numerical_damping_factor_;
    DiscreteVariable<Vecd> *dv_force_prior_;
};

template <typename...>
class StructureIntegration2ndHalf;

template <typename... Parameters>
class StructureIntegration2ndHalf<Inner<OneLevel, Parameters...>>
    : public Interaction<Inner<Parameters...>>, public StructureDynamicsVariables
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    explicit StructureIntegration2ndHalf(Inner<Parameters...> &inner_relation);
    virtual ~StructureIntegration2ndHalf() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *pos_, *vel_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *Vol0_;
        Vecd *vel_;
        Matd *B_, *dF_dt_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Matd *F_, *dF_dt_;
    };
};

template <typename...>
class StructureNumericalDamping;

template <class MaterialType, typename... Parameters>
class StructureNumericalDamping<Inner<WithUpdate, MaterialType, Parameters...>>
    : public Interaction<Inner<Parameters...>>,
      public StructureDynamicsVariables,
      public ForcePriorCK
{
    using BaseInteraction = Interaction<Inner<Parameters...>>;
    using Adaptation = typename Inner<Parameters...>::SourceType::Adaptation;
    using SmoothingLengthRatioType = typename Adaptation::SmoothingLengthRatioType;
    using ConstituteKernel = typename MaterialType::ConstituteKernel;

  public:
    explicit StructureNumericalDamping(Inner<Parameters...> &inner_relation, Real numerical_damping_factor = 0.125);
    virtual ~StructureNumericalDamping() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        ConstituteKernel constitute_;
        SmoothingLengthRatioType h_ratio_;
        Vecd zero_;
        Real h_ref_, numerical_damping_factor_;
        Real *Vol0_;
        Vecd *pos_, *vel_, *numerical_damping_force_;
        Matd *F_;
    };

  protected:
    MaterialType &material_;
    Adaptation &adaptation_;
    Real h_ref_, numerical_damping_factor_;
    DiscreteVariable<Vecd> *dv_numerical_damping_force_;
};

} // namespace solid_dynamics
} // namespace SPH
#endif // STRUCTURE_DYNAMICS_H
