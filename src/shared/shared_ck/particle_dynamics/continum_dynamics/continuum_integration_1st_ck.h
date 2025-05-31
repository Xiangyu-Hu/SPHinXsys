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
 * @file 	continuum_integration_1st_ck.h
 * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
 * @details CK and SYCL version.
 * @author	Shuang Li,Xiangyu Hu and Shuaihao Zhang
 */
#ifndef CONTINUUM_INTEGRATION_1ST_CK_H
#define CONTINUUM_INTEGRATION_1ST_CK_H

#include "acoustic_step_1st_half.h"
#include "constraint_dynamics.h"
#include "general_continuum.h"
#include "general_continuum.hpp"
namespace SPH
{
namespace continuum_dynamics
{

template <class BaseInteractionType>
class PlasticAcousticStep : public fluid_dynamics::AcousticStep<BaseInteractionType>
{

  public:
    template <class DynamicsIdentifier>
    explicit PlasticAcousticStep(DynamicsIdentifier &identifier);
    virtual ~PlasticAcousticStep() {};

  protected:
    PlasticContinuum &plastic_continuum_;
    // DiscreteVariable<Vecd> *dv_pos_;
    DiscreteVariable<Mat3d> *dv_stress_tensor_3D_, *dv_strain_tensor_3D_, *dv_stress_rate_3D_, *dv_strain_rate_3D_;
    DiscreteVariable<Matd> *dv_velocity_gradient_;
};

template <typename...>
class PlasticAcousticStep1stHalf;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStep<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = PlasticAcousticStep<Interaction<Inner<Parameters...>>>;

  public:
    explicit PlasticAcousticStep1stHalf(Inner<Parameters...> &inner_relation);
    virtual ~PlasticAcousticStep1stHalf() {};

    class InitializeKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void initialize(size_t index_i, Real dt = 0.0);

      protected:
        Real *rho_, *p_, *drho_dt_;
        Vecd *vel_, *dpos_;
        Mat3d *stress_tensor_3D_;
    };

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        KernelCorrectionType correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *p_, *drho_dt_, *mass_;
        Vecd *force_;
        Mat3d *stress_tensor_3D_;
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
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStep1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStep<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using BaseInteraction = PlasticAcousticStep<Interaction<Contact<Parameters...>>>;

  public:
    explicit PlasticAcousticStep1stHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~PlasticAcousticStep1stHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        KernelCorrectionType correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *mass_, *p_, *drho_dt_;
        Vecd *force_, *force_prior_;
        Mat3d *stress_tensor_3D_;

        Real *wall_Vol_;
        Vecd *wall_acc_ave_;
    };

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

using PlasticAcousticStep1stHalfWithWallRiemannCK =
    PlasticAcousticStep1stHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrection>,
                               Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrection>>;
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_INTEGRATION_1ST_CK_H