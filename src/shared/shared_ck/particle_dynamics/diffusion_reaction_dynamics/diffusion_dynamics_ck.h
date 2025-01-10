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
 * @file    diffusion_dynamics_ck.h
 * @brief   tbd.
 * @author  Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_CK_H
#define DIFFUSION_DYNAMICS_CK_H

#include "diffusion_reaction.h"
#include "interaction_algorithms_ck.hpp"
#include "interaction_ck.hpp"
#include "sphinxsys_variable_array.h"

namespace SPH
{
class CorrectedKernelGradientInnerCK
{
    DiscreteVariable<Matd> *dv_B_;

  public:
    explicit CorrectedKernelGradientInnerCK(BaseParticles *particles)
        : dv_B_(particles->getVariableByName<Matd>("LinearGradientCorrectionMatrix")) {};

    class ComputingKernel
    {
        Matd *B_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : B_(encloser.dv_B_->DelegatedData(ex_policy)){};

        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return 0.5 * dW_ijV_j * (B_[index_i] + B_[index_j]) * e_ij;
        };
    };
};

template <typename... InteractionTypes>
class DiffusionRelaxationCK;

template <class DiffusionType, class BaseInteractionType>
class DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>
    : public BaseInteractionType
{
  protected:
    UniquePtrsKeeper<Entity> variable_array_ptrs_;

  public:
    template <class DynamicsIdentifier>
    explicit DiffusionRelaxationCK(DynamicsIdentifier &identifier);
    virtual ~DiffusionRelaxationCK() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_dt_;
        UnsignedInt number_of_species_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_, *diffusion_dt_;
        UnsignedInt number_of_species_;
    };

  protected:
    StdVec<DiffusionType *> diffusions_;
    DiscreteVariableArray<Real> *dv_diffusion_species_array_;
    DiscreteVariableArray<Real> *dv_gradient_species_array_;
    DiscreteVariableArray<Real> *dv_diffusion_dt_array_;
    Real smoothing_length_sq_;

  private:
    void getDiffusions();
};

template <class DiffusionType, class BaseInteractionType>
class DiffusionRelaxationCK<RungeKutta, DiffusionType, BaseInteractionType>
    : public DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>
{
    using BaseDynamicsType = DiffusionRelaxationCK<ForwardEuler, DiffusionType, BaseInteractionType>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

  protected:
    DiscreteVariableArray<Real> *dv_diffusion_species_array_s_;
};

template <class DiffusionType, class BaseInteractionType>
class DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType, BaseInteractionType>
    : public DiffusionRelaxationCK<RungeKutta, DiffusionType, BaseInteractionType>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RungeKutta, DiffusionType, BaseInteractionType>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InitializeKernel : public BaseDynamicsType::InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_s_;
        DataArray<Real> *diffusion_species_;
    };
};

template <class DiffusionType, class BaseInteractionType>
class DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType, BaseInteractionType>
    : public DiffusionRelaxationCK<RungeKutta, DiffusionType, BaseInteractionType>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RungeKutta, DiffusionType, BaseInteractionType>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class UpdateKernel : public BaseDynamicsType::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_s_;
    };
};

template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
class DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>
    : public DiffusionRelaxationCK<TimeSteppingType, DiffusionType, Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = DiffusionRelaxationCK<TimeSteppingType, DiffusionType, Interaction<Inner<Parameters...>>>;
    using GradientKernel = typename KernelGradientType::ComputingKernel;
    using InterParticleDiffusionCoeff = typename DiffusionType::InterParticleDiffusionCoeff;

  public:
    explicit DiffusionRelaxationCK(Relation<Inner<Parameters...>> &inner_relation);
    virtual ~DiffusionRelaxationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        InterParticleDiffusionCoeff *inter_particle_diffusion_coeff_;
        DataArray<Real> *diffusion_species_, *gradient_species_, *diffusion_dt_;
        UnsignedInt number_of_species_;
        GradientKernel gradient_;
        Real *Vol_;
        Real smoothing_length_sq_;
    };

  protected:
    KernelGradientType kernel_gradient_;
    ConstantArray<DiffusionType, InterParticleDiffusionCoeff> ca_inter_particle_diffusion_coeff_;
    DiscreteVariable<Real> *dv_Vol_;
};
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_H
