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
 * @file 	acoustic_step_1st_half.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef ACOUSTIC_STEP_1ST_HALF_H
#define ACOUSTIC_STEP_1ST_HALF_H

#include "base_fluid_dynamics.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"
#include "riemann_solver_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{

template <class BaseInteractionType>
class AcousticStep : public BaseInteractionType
{

  public:
    template <class DynamicsIdentifier>
    explicit AcousticStep(DynamicsIdentifier &identifier);
    virtual ~AcousticStep() {};

  protected:
    DiscreteVariable<Real> *dv_Vol_, *dv_rho_, *dv_mass_, *dv_p_, *dv_drho_dt_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_dpos_, *dv_force_, *dv_force_prior_;
};

template <typename...>
class AcousticStep1stHalf;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Inner<Parameters...>>>
{
    using FluidType = typename RiemannSolverType::SourceFluid;
    using EosKernel = typename FluidType::EosKernel;
    using BaseInteraction = AcousticStep<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep1stHalf(Inner<Parameters...> &inner_relation);
    virtual ~AcousticStep1stHalf() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        EosKernel eos_;
        Real *rho_, *p_, *drho_dt_;
        Vecd *vel_, *dpos_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *p_, *drho_dt_;
        Vecd *force_;
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
    KernelCorrectionType kernel_correction_;
    FluidType &fluid_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using FluidType = typename RiemannSolverType::SourceFluid;
    using BaseInteraction = AcousticStep<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep1stHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~AcousticStep1stHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *mass_, *p_, *drho_dt_;
        Vecd *vel_, *force_, *force_prior_;
        Real *wall_Vol_;
        Vecd *wall_acc_ave_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    FluidType &fluid_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Contact<Parameters...>>>
{
    using SourceFluidType = typename RiemannSolverType::SourceFluid;
    using TargetFluidType = typename RiemannSolverType::TargetFluid;
    using BaseInteraction = AcousticStep<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep1stHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~AcousticStep1stHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        CorrectionKernel contact_correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *p_, *drho_dt_;
        Vecd *force_;
        Real *contact_Vol_, *contact_p_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    StdVec<KernelCorrectionType> contact_kernel_corrections_;
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_, dv_contact_p_;
};

using AcousticStep1stHalfWithWallRiemannCK =
    AcousticStep1stHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>,
                        Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>;
using AcousticStep1stHalfWithWallRiemannCorrectionCK =
    AcousticStep1stHalf<Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>,
                        Contact<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_1ST_HALF_H
