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
class KernelGradientInnerCK
{
  public:
    explicit KernelGradientInnerCK(BaseParticles *inner_particles) {};
    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser){};
        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return dW_ijV_j * e_ij;
        };
    };
};

class CorrectedKernelGradientInnerCK
{
    DiscreteVariable<Matd> *dv_B_;

  public:
    explicit CorrectedKernelGradientInnerCK(BaseParticles *particles)
        : dv_B_(particles->getVariableByName<Matd>("LinearCorrectionMatrix")) {};

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

class KernelGradientContactCK
{
  public:
    KernelGradientContactCK(BaseParticles *inner_particles, BaseParticles *contact_particles) {};
    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser){};
        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return dW_ijV_j * e_ij;
        };
    };
};

class CorrectedKernelGradientContactCK
{
    DiscreteVariable<Matd> *dv_B_;
    DiscreteVariable<Matd> *dv_contact_B_;

  public:
    CorrectedKernelGradientContactCK(BaseParticles *inner_particles, BaseParticles *contact_particles)
        : dv_B_(inner_particles->getVariableByName<Matd>("LinearCorrectionMatrix")),
          dv_contact_B_(contact_particles->getVariableByName<Matd>("LinearCorrectionMatrix")) {};

    class ComputingKernel
    {
        Matd *B_;
        Matd *contact_B_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : B_(encloser.dv_B_->DelegatedData(ex_policy)),
              contact_B_(encloser.dv_contact_B_->DelegatedData(ex_policy)){};

        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return 0.5 * dW_ijV_j * (B_[index_i] + contact_B_[index_j]) * e_ij;
        };
    };
};

template <typename... InteractionTypes>
class DiffusionRelaxationCK;

template <class DiffusionType>
class DiffusionRelaxationCK<ForwardEuler, DiffusionType>
{
    StdVec<DiffusionType *> getConcreteDiffusions(AbstractDiffusion &abstract_diffusion);
    StdVec<DiscreteVariable<Real> *> getDiffusionSpecies(BaseParticles &particles);
    StdVec<DiscreteVariable<Real> *> getGradientSpecies(BaseParticles &particles);
    StdVec<DiscreteVariable<Real> *> getSpeciesChangeRates(BaseParticles &particles);

  public:
    DiffusionRelaxationCK(SPHBody &sph_body, AbstractDiffusion &abstract_diffusion);
    explicit DiffusionRelaxationCK(SPHBody &sph_body);
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
    DiscreteVariableArray<Real> dv_diffusion_species_array_;
    DiscreteVariableArray<Real> dv_gradient_species_array_;
    DiscreteVariableArray<Real> dv_diffusion_dt_array_;
    Real smoothing_length_sq_;
};

template <class DiffusionType>
class DiffusionRelaxationCK<RungeKutta, DiffusionType>
    : public DiffusionRelaxationCK<ForwardEuler, DiffusionType>
{
    using BaseDynamicsType = DiffusionRelaxationCK<ForwardEuler, DiffusionType>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

  protected:
    DiscreteVariableArray<Real> dv_diffusion_species_array_s_;

  private:
    StdVec<DiscreteVariable<Real> *> getIntermediateSpecies(BaseParticles &particles);
};

template <class DiffusionType>
class DiffusionRelaxationCK<RungeKutta1stStage, DiffusionType>
    : public DiffusionRelaxationCK<RungeKutta, DiffusionType>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RungeKutta, DiffusionType>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InitializeKernel : public BaseDynamicsType::InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_;
        DataArray<Real> *diffusion_species_s_;
    };
};

template <class DiffusionType>
class DiffusionRelaxationCK<RungeKutta2ndStage, DiffusionType>
    : public DiffusionRelaxationCK<RungeKutta, DiffusionType>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RungeKutta, DiffusionType>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args);
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

template <class TimeSteppingType, class DiffusionType, class KernelGradientType>
class DiffusionRelaxationCK<TimeSteppingType, DiffusionType, KernelGradientType>
    : public DiffusionRelaxationCK<TimeSteppingType, DiffusionType>
{
    using BaseRelaxation = DiffusionRelaxationCK<TimeSteppingType, DiffusionType>;
    using GradientKernel = typename KernelGradientType::ComputingKernel;
    using InterParticleDiffusionCoeff = typename DiffusionType::InterParticleDiffusionCoeff;

  public:
    template <typename... Args>
    explicit DiffusionRelaxationCK(SPHBody &sph_body, Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

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

template <class TimeSteppingType, class DiffusionType, class KernelGradientType, typename... Parameters>
class DiffusionRelaxationCK<Inner<OneLevel, TimeSteppingType, DiffusionType, KernelGradientType, Parameters...>>
    : public DiffusionRelaxationCK<TimeSteppingType, DiffusionType, KernelGradientType>,
      public Interaction<Inner<Parameters...>>
{
    using BaseRelaxation = DiffusionRelaxationCK<TimeSteppingType, DiffusionType, KernelGradientType>;
    using BaseInteraction = Interaction<Inner<Parameters...>>;

  public:
    template <typename... Args>
    explicit DiffusionRelaxationCK(Relation<Inner<Parameters...>> &relation, Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InteractKernel : public BaseRelaxation::InteractKernel, public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);
    };
};
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_H
