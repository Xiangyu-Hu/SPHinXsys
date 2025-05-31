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
 * @file 	continuum_integration.h
 * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
 * @details CK and SYCL version.
 * @author	Shuang Li, Xiangyu Hu and Shuaihao Zhang
 */
#ifndef CONTINUUM_INTEGRATION_2ND_CK_H
#define CONTINUUM_INTEGRATION_2ND_CK_H

#include "continuum_integration_1st_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
// step2 inner
template <typename...>
class PlasticAcousticStep2ndHalf;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStep2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStep<Interaction<Inner<Parameters...>>>
{
    using PlasticKernel = typename PlasticContinuum::PlasticKernel;
    using BaseInteraction = PlasticAcousticStep<Interaction<Inner<Parameters...>>>;

  public:
    explicit PlasticAcousticStep2ndHalf(Inner<Parameters...> &inner_relation);
    virtual ~PlasticAcousticStep2ndHalf() {};

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
        KernelCorrectionType correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *drho_dt_;
        Vecd *vel_, *force_;
        Matd *velocity_gradient_;
    };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        PlasticKernel plastic_kernel_;
        Real *rho_, *drho_dt_;
        Matd *velocity_gradient_;
        Mat3d *stress_tensor_3D_, *strain_tensor_3D_, *stress_rate_3D_, *strain_rate_3D_;
    };

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStep2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStep<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using BaseInteraction = PlasticAcousticStep<Interaction<Contact<Parameters...>>>;

  public:
    explicit PlasticAcousticStep2ndHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~PlasticAcousticStep2ndHalf() {};

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        KernelCorrectionType correction_;
        RiemannSolverType riemann_solver_;
        Real *Vol_, *rho_, *drho_dt_;
        Vecd *vel_, *force_;
        Real *wall_Vol_;
        Vecd *wall_vel_ave_, *wall_n_;
        Matd *velocity_gradient_;
    };

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

using PlasticAcousticStep2ndHalfWithWallRiemannCK =
    PlasticAcousticStep2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrection>,
                               Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrection>>;
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_INTEGRATION_2ND_CK_H