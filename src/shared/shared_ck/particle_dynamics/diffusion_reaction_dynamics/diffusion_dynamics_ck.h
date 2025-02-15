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
#include "interaction_ck.hpp"
#include "sphinxsys_variable_array.h"
#include "sphinxsys_constant.h"

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
        Vecd operator()(UnsignedInt index_i, UnsignedInt index_j, Real dW_ijV_j, const Vecd &e_ij)
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

        Vecd operator()(UnsignedInt index_i, UnsignedInt index_j, Real dW_ijV_j, const Vecd &e_ij)
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
        Vecd operator()(UnsignedInt index_i, UnsignedInt index_j, Real dW_ijV_j, const Vecd &e_ij)
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

        Vecd operator()(UnsignedInt index_i, UnsignedInt index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return 0.5 * dW_ijV_j * (B_[index_i] + contact_B_[index_j]) * e_ij;
        };
    };
};

template <typename... InteractionTypes>
class DiffusionRelaxationCK;

template <class DiffusionType, class BaseInteractionType>
class DiffusionRelaxationCK<DiffusionType, BaseInteractionType>
    : public BaseInteractionType
{
    StdVec<DiffusionType *> getConcreteDiffusions(AbstractDiffusion &abstract_diffusion);

  public:
    template <class DynamicsIdentifier>
    DiffusionRelaxationCK(DynamicsIdentifier &identifier, AbstractDiffusion *abstract_diffusion);
    template <class DynamicsIdentifier>
    explicit DiffusionRelaxationCK(DynamicsIdentifier &identifier);
    template <typename BodyRelationType, typename FirstArg>
    DiffusionRelaxationCK(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : DiffusionRelaxationCK(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~DiffusionRelaxationCK() {};
    StdVec<DiffusionType *> &getDiffusions() { return diffusions_; };
    DiscreteVariableArray<Real> &dvGradientSpeciesArray() { return dv_gradient_species_array_; };
    StdVec<DiscreteVariable<Real> *> getDiffusionVariables(BaseParticles *particles, const std::string &suffix);
    StdVec<DiscreteVariable<Real> *> getGradientVariables(BaseParticles *particles, const std::string &suffix);

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        DataArray<Real> *diffusion_species_, *gradient_species_, *diffusion_dt_;
        UnsignedInt number_of_species_;
    };

  protected:
    StdVec<DiffusionType *> diffusions_;
    DiscreteVariableArray<Real> dv_diffusion_species_array_;
    DiscreteVariableArray<Real> dv_gradient_species_array_;
    DiscreteVariableArray<Real> dv_diffusion_dt_array_;
};

template <class DiffusionType, class KernelGradientType, class... Parameters>
class DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelGradientType, Parameters...>>
    : public DiffusionRelaxationCK<DiffusionType, Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = DiffusionRelaxationCK<DiffusionType, Interaction<Inner<Parameters...>>>;
    using GradientKernel = typename KernelGradientType::ComputingKernel;
    using InterParticleDiffusionCoeff = typename DiffusionType::InterParticleDiffusionCoeff;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        GradientKernel gradient_;
        InterParticleDiffusionCoeff *inter_particle_diffusion_coeff_;
        Real *Vol_;
        Real smoothing_length_sq_;
    };

  protected:
    KernelGradientType kernel_gradient_;
    ConstantArray<DiffusionType, InterParticleDiffusionCoeff> ca_inter_particle_diffusion_coeff_;
    DiscreteVariable<Real> *dv_Vol_;
    Real smoothing_length_sq_;
};

template <class DiffusionType, template <typename> class BoundaryType, class KernelGradientType>
class DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType>, KernelGradientType>>
    : public DiffusionRelaxationCK<DiffusionType, Interaction<Contact<>>>
{
    UniquePtrsKeeper<DiscreteVariableArray<Real>> contact_transfer_array_ptrs_keeper_;
    UniquePtrsKeeper<KernelGradientType> kernel_gradient_ptrs_keeper_;
    UniquePtrsKeeper<BoundaryType<DiffusionType>> boundary_ptrs_keeper_;
    using BaseInteraction = DiffusionRelaxationCK<DiffusionType, Interaction<Contact<>>>;
    using GradientKernel = typename KernelGradientType::ComputingKernel;
    using BoundaryKernel = typename BoundaryType<DiffusionType>::ComputingKernel;

  public:
    template <typename... Args>
    explicit DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser,
                       UnsignedInt contact_index);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        Real *contact_Vol_;
        DataArray<Real> *contact_transfer_;
        GradientKernel gradient_;
        BoundaryKernel boundary_flux_;
    };

  protected:
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariableArray<Real> *> contact_dv_transfer_array_;
    StdVec<KernelGradientType *> contact_kernel_gradient_method_;
    StdVec<BoundaryType<DiffusionType> *> contact_boundary_method_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
    : public DiffusionRelaxationCK<RelationType<InteractionOnly, InteractionParameters...>>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RelationType<InteractionOnly, InteractionParameters...>>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args) : BaseDynamicsType(std::forward<Args>(args)...){};
    virtual ~DiffusionRelaxationCK() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_dt_;
        UnsignedInt number_of_species_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_, *diffusion_dt_;
        UnsignedInt number_of_species_;
    };
};

template <template <typename...> class RelationType, class... InteractionParameters>
class DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta1stStage, InteractionParameters...>>
    : public DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InitializeKernel : public BaseDynamicsType::InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_;
        DataArray<Real> *diffusion_species_s_;
    };

  protected:
    DiscreteVariableArray<Real> dv_diffusion_species_array_s_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class DiffusionRelaxationCK<RelationType<OneLevel, RungeKutta2ndStage, InteractionParameters...>>
    : public DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class UpdateKernel : public BaseDynamicsType::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0);

      protected:
        DataArray<Real> *diffusion_species_s_;
    };

  protected:
    DiscreteVariableArray<Real> dv_diffusion_species_array_s_;
};

template <class DiffusionType>
class Dirichlet<DiffusionType>
{
    using InterParticleDiffusionCoeff = typename DiffusionType::InterParticleDiffusionCoeff;

  public:
    template <class DiffusionDynamics>
    Dirichlet(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Vecd operator()(UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j, const Vecd &e_ij, const Vecd &vec_r_ij);

      protected:
        Real smoothing_length_sq_;
        DataArray<Real> *gradient_species_;
        DataArray<Real> *contact_gradient_species_;
        InterParticleDiffusionCoeff *inter_particle_diffusion_coeff_;
    };

  protected:
    Real smoothing_length_sq_;
    DiscreteVariableArray<Real> &dv_gradient_species_array_;
    DiscreteVariableArray<Real> contact_dv_gradient_species_array_;
    ConstantArray<DiffusionType, InterParticleDiffusionCoeff> ca_inter_particle_diffusion_coeff_;
};

template <class DiffusionType>
class Neumann<DiffusionType>
{
  public:
    template <class DiffusionDynamics>
    Neumann(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Vecd operator()(UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j, const Vecd &e_ij, const Vecd &vec_r_ij);

      protected:
        DataArray<Real> *contact_species_flux_;
    };

  protected:
    DiscreteVariableArray<Real> contact_dv_species_flux_array_;
};
} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_H
