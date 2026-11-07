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
    virtual ~AcousticStep(){};

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_p_, *dv_drho_dt_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_dpos_, *dv_force_, *dv_force_prior_;
    /**
     * Optional background hydrostatic pressure p_h (Doc87 R3-A / Doc97–98 / Doc105).
     * Always registered (default 0). Initialize sets p = EOS(ρ) + BackgroundPressure.
     *
     * Momentum pressure is selected by PerturbationMomentum ∈ {0,1} (default 0):
     *   0 → total p (keep GravityForceCK to cancel discrete ∇p_h when using p_h)
     *   1 → p' only (closed-box; omit GravityForceCK, keep buoyancy)
     *
     * Density dissipation mode AcousticDensMode ∈ {0,1,2,3} (default 0 = legacy):
     *   0 → DissipativeUJump on the same Adami pressure used for momentum
     *       (wall: p_i_mom − (p_i_mom+C_ext); restores non-zero wall dens when C≠0).
     *       Inner uses total p. This is the pre-Doc98 wall dens semantics.
     *   1 → DissipativeUJump on p' = p − BackgroundPressure. Wall reference ghost
     *       increment uses ReferenceBodyAcceleration (independent of ForcePrior),
     *       so dens jump = C_ext − C_ref (zero only when C_ext==C_ref).
     *       Closed-box P2b (g absorbed in BackgroundPressure, GravityForceCK omitted):
     *       set ReferenceBodyAcceleration=0 so C_ref=0; otherwise C_ext=0 & C_ref=C_g
     *       injects spurious wall dens. Set C_ref from g only when g remains in ForcePrior.
     *   2 → Doc186: dens as mode 1, but wall *momentum* Adami omits C_ext (p_j_mom=p_i_mom).
     *       Use when body-force hydrostatic already lives in p'/ρ (e.g. uniform-buoyancy
     *       balance init): mode 1 then double-counts ForcePrior at walls and sustains drift.
     *       Keep --cref-buoyancy so dens C_ext≈C_ref. Not a default for free convection
     *       without consistent p' hydro.
     *   3 → Doc201: momentum as mode 1; *omit* pressure DissipativeUJump dens (inner+wall).
     *       Continuity dens then comes only from 2nd-half velocity divergence.
     *       Needed when p' itself carries stratified hydrostatic (Doc200 strat_init): mode 1
     *       DissipativeUJump(p'_i−p'_j) is nonzero at rest along ∇p' and seeds drift.
     *       Mode 2 is not a substitute (omitting wall C_ext made Doc201 B1 residual worse).
     */
    DiscreteVariable<Real> *dv_p_background_;
    DiscreteVariable<Real> *dv_perturbation_momentum_;
    DiscreteVariable<Real> *dv_acoustic_dens_mode_;
    DiscreteVariable<Vecd> *dv_reference_body_acceleration_;
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
    using CorrectionDataType = typename KernelCorrectionType::CorrectionDataType;

  public:
    explicit AcousticStep1stHalf(Inner<Parameters...> &inner_relation);
    virtual ~AcousticStep1stHalf(){};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        EosKernel eos_;
        Real *rho_, *p_, *drho_dt_, *p_background_;
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
        Real *Vol_, *rho_, *p_, *p_background_, *perturbation_momentum_, *acoustic_dens_mode_, *drho_dt_;
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
    virtual ~AcousticStep1stHalf(){};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        CorrectionKernel correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *mass_, *p_, *p_background_, *perturbation_momentum_, *acoustic_dens_mode_,
            *drho_dt_;
        Vecd *vel_, *force_, *force_prior_;
        Vecd *reference_body_acceleration_;
        Real *contact_Vol_;
        Vecd *wall_acc_ave_, *wall_n_;
        Real *wall_free_slip_; /**< Doc201: >0.5 → project mom force onto wall normal */
        Real *wall_contact_omit_; /**< Doc202: >0.5 → skip both acoustic halves */
        Real *wall_contact_omit_1st_; /**< Doc202 H1c: >0.5 → skip 1st-half pair only */
    };

  protected:
    KernelCorrectionType kernel_correction_;
    FluidType &fluid_;
    RiemannSolverType riemann_solver_;
    /** Shared with AcousticStep2ndHalf Contact<Wall>; default 0 = legacy. */
    StdVec<DiscreteVariable<Real> *> dv_wall_free_slip_;
    /** Doc202 H1: omit both halves; shared register with 2nd-half. */
    StdVec<DiscreteVariable<Real> *> dv_wall_contact_omit_;
    /** Doc202 H1c: omit 1st-half only; shared register. */
    StdVec<DiscreteVariable<Real> *> dv_wall_contact_omit_1st_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class AcousticStep1stHalf<Contact<RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public AcousticStep<Interaction<Contact<Parameters...>>>
{
    using SourceFluidType = typename RiemannSolverType::SourceFluid;
    using TargetFluidType = typename RiemannSolverType::TargetFluid;
    using BaseInteraction = AcousticStep<Interaction<Contact<Parameters...>>>;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using CorrectionDataType = typename KernelCorrectionType::CorrectionDataType;

  public:
    explicit AcousticStep1stHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~AcousticStep1stHalf(){};

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
    StdVec<DiscreteVariable<Real> *> dv_contact_p_;
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
