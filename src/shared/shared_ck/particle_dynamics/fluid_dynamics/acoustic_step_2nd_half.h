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
 * @file 	acoustic_step_2nd_half.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef ACOUSTIC_STEP_2ND_HALF_H
#define ACOUSTIC_STEP_2ND_HALF_H

#include "acoustic_step_1st_half.h"

namespace SPH
{
namespace fluid_dynamics
{

template <typename...>
class AcousticStep2ndHalf;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Inner<Parameters...>>>
{
    using FluidType = typename RiemannSolverType::SourceFluid;
    using BaseInteraction = AcousticStep<Interaction<Inner<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep2ndHalf(Inner<Parameters...> &inner_relation);
    virtual ~AcousticStep2ndHalf() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
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
        Real *Vol_, *rho_, *drho_dt_;
        Vecd *vel_, *force_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real *rho_, *drho_dt_;
    };

  protected:
    template <typename...>
    friend class FSI::PressureForceFromFluid;

    KernelCorrectionType kernel_correction_;
    FluidType &fluid_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using FluidType = typename RiemannSolverType::SourceFluid;
    using BaseInteraction = AcousticStep<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep2ndHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~AcousticStep2ndHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *drho_dt_;
        Vecd *vel_, *force_;
        Real *wall_Vol_;
        Vecd *wall_vel_ave_, *wall_n_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    FluidType &fluid_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep2ndHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Contact<Parameters...>>>
{
    using SourceFluidType = typename RiemannSolverType::SourceFluid;
    using TargetFluidType = typename RiemannSolverType::TargetFluid;
    using BaseInteraction = AcousticStep<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit AcousticStep2ndHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~AcousticStep2ndHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *drho_dt_;
        Vecd *vel_, *force_;
        Real *contact_Vol_;
        Vecd *contact_vel_;
    };

  protected:
    KernelCorrectionType kernel_correction_;
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<DiscreteVariable<Real> *> dv_contact_Vol_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_vel_;
};

using AcousticStep2ndHalfWithWallNoRiemannCK =
    AcousticStep2ndHalf<Inner<OneLevel, NoRiemannSolverCK, NoKernelCorrectionCK>,
                        Contact<Wall, NoRiemannSolverCK, NoKernelCorrectionCK>>;
using AcousticStep2ndHalfWithWallRiemannCK =
    AcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrectionCK>,
                        Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrectionCK>>;
using AcousticStep2ndHalfWithWallRiemannCorrectionCK =
    AcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, LinearCorrectionCK>,
                        Contact<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // ACOUSTIC_STEP_2ND_HALF_H
