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
 * @file    diffusion_dynamics_ck.h
 * @brief   tbd.
 * @author  Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_CK_H
#define DIFFUSION_DYNAMICS_CK_H

#include "diffusion_reaction.h"
#include "interaction_ck.hpp"
#include "sphinxsys_constant.h"

#include <string>
#include <tuple>
#include <utility>

namespace SPH
{
template <class DiffusionType>
StdVec<DiffusionType *> obtainConcreteDiffusions(AbstractDiffusion &abstract_diffusion);

template <class DiffusionType>
StdVec<std::string> obtainSpeciesNames(StdVec<DiffusionType *> &diffusions);

template <class DiffusionType, class BaseLocalDynamicsType>
class BaseDiffusionDynamics : public BaseLocalDynamicsType
{
  public:
    typedef DiffusionType Diffusion;
    template <class DynamicsIdentifier>
    BaseDiffusionDynamics(DynamicsIdentifier &identifier, AbstractDiffusion *abstract_diffusion);
    template <class DynamicsIdentifier>
    explicit BaseDiffusionDynamics(DynamicsIdentifier &identifier);
    template <typename BodyRelationType, typename FirstArg>
    BaseDiffusionDynamics(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : BaseDiffusionDynamics(parameters.identifier_, std::get<0>(parameters.others_)) {};
    virtual ~BaseDiffusionDynamics() {};
    StdVec<DiffusionType *> &getDiffusions() { return diffusions_; };
    StdVec<std::string> &getSpeciesNames() { return species_names_; };
    DiscreteVariable<Real> *registerSpecies(
        BaseParticles *particles, StdVec<std::string> &species_names, std::string suffix = "");
    DiscreteVariable<Real> *getSpeciesByName(
        BaseParticles *particles, StdVec<std::string> &species_names, std::string suffix = "");

  protected:
    StdVec<DiffusionType *> diffusions_;
    StdVec<std::string> species_names_;
    DiscreteVariable<Real> *dv_species_;
    DiscreteVariable<Real> *dv_species_dt_;
};

template <typename... InteractionTypes>
class DiffusionRelaxationCK;

template <class DiffusionType, class KernelCorrectionType, class... Parameters>
class DiffusionRelaxationCK<Inner<InteractionOnly, DiffusionType, KernelCorrectionType, Parameters...>>
    : public BaseDiffusionDynamics<DiffusionType, Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = BaseDiffusionDynamics<DiffusionType, Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using InterParticleCoeff = typename DiffusionType::InterParticleCoeff;

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
        CorrectionKernel correction_;
        MultiEntryView<Real> species_, species_dt_;
        InterParticleCoeff *inter_particle_diffusion_coeff_;
        Real *Vol_;
        Real smoothing_length_sq_;
    };

  protected:
    KernelCorrectionType kernel_correction_method_;
    ComputingKernelArray<DiffusionType, InterParticleCoeff> cka_inter_particle_diffusion_coeff_;
    Real smoothing_length_sq_;
};

template <class DiffusionType, class... OtherTypes, template <typename...> class BoundaryType,
          class KernelCorrectionType, class... Parameters>
class DiffusionRelaxationCK<Contact<InteractionOnly, BoundaryType<DiffusionType, OtherTypes...>,
                                    KernelCorrectionType, Parameters...>>
    : public BaseDiffusionDynamics<DiffusionType, Interaction<Contact<Parameters...>>>
{
    UniquePtrKeeper<BoundaryType<DiffusionType, OtherTypes...>> boundary_keeper_;
    using BaseInteraction = BaseDiffusionDynamics<DiffusionType, Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using BoundaryKernel = typename BoundaryType<DiffusionType, OtherTypes...>::ComputingKernel;

  public:
    template <typename... Args>
    explicit DiffusionRelaxationCK(
        Contact<Parameters...> &contact_relation, AbstractDiffusion *abstract_diffusion, Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(UnsignedInt index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        MultiEntryView<Real> species_, species_dt_;
        Real *contact_Vol_;
        MultiEntryView<Real> contact_transfer_;
        BoundaryKernel boundary_flux_;
    };

  protected:
    KernelCorrectionType kernel_correction_method_;
    DiscreteVariable<Real> *contact_dv_transfer_;
    BoundaryType<DiffusionType, OtherTypes...> *contact_boundary_method_;
};

template <template <typename...> class RelationType, class... InteractionParameters>
class DiffusionRelaxationCK<RelationType<OneLevel, ForwardEuler, InteractionParameters...>>
    : public DiffusionRelaxationCK<RelationType<InteractionOnly, InteractionParameters...>>
{
    using BaseDynamicsType = DiffusionRelaxationCK<RelationType<InteractionOnly, InteractionParameters...>>;
    using DiffusionType = typename BaseDynamicsType::Diffusion;
    using InverseVolumetricCapacity = typename DiffusionType::InverseVolumetricCapacity;

  public:
    template <typename... Args>
    DiffusionRelaxationCK(Args &&...args);
    virtual ~DiffusionRelaxationCK() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(UnsignedInt index_i, Real dt = 0.0);

      protected:
        MultiEntryView<Real> species_dt_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(UnsignedInt index_i, Real dt = 0.0);

      protected:
        MultiEntryView<Real> species_, species_dt_;
        InverseVolumetricCapacity *cv1_;
    };

  protected:
    ComputingKernelArray<DiffusionType, InverseVolumetricCapacity> cka_inverse_volume_capacity_;
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
        MultiEntryView<Real> species_, species_s_;
    };

  protected:
    DiscreteVariable<Real> *dv_species_s_;
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
        MultiEntryView<Real> species_, species_s_;
    };

  protected:
    DiscreteVariable<Real> *dv_species_s_;
};

template <class DiffusionType>
class Dirichlet<DiffusionType>
{
    using InterParticleCoeff = typename DiffusionType::InterParticleCoeff;

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
        MultiEntryView<Real> species_, contact_species_;
        InterParticleCoeff *inter_particle_diffusion_coeff_;
    };

  protected:
    Real smoothing_length_sq_;
    DiscreteVariable<Real> *dv_species_;
    DiscreteVariable<Real> *dv_contact_species_;
    ComputingKernelArray<DiffusionType, InterParticleCoeff> cka_inter_particle_diffusion_coeff_;
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
        Vecd *contact_n_;
        MultiEntryView<Real> contact_species_flux_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_contact_n_;
    DiscreteVariable<Real> *dv_contact_species_flux_;
};

template <class SourceDiffusionType, class TargetDiffusionType>
class Robin<SourceDiffusionType, TargetDiffusionType>
{
  public:
    template <class DiffusionDynamics>
    Robin(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles,
          AbstractDiffusion *abstract_diffusion, Real transfer_coefficient);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Vecd operator()(UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j, const Vecd &e_ij, const Vecd &vec_r_ij);

      protected:
        Real smoothing_length_sq_;
        MultiEntryView<Real> species_, contact_species_;
        Vecd *contact_n_;
        Real transfer_coefficient_;
    };

  protected:
    StdVec<TargetDiffusionType *> target_diffusions_;
    Real smoothing_length_sq_;
    DiscreteVariable<Real> *dv_species_;
    StdVec<std::string> target_species_names_;
    DiscreteVariable<Real> *dv_contact_species_;
    DiscreteVariable<Vecd> *dv_contact_n_;
    Real transfer_coefficient_;
};

} // namespace SPH
#endif // DIFFUSION_DYNAMICS_CK_H
