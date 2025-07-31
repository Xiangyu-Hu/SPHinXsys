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
 * @file 	erosion_dynamics_2nd_ck.h
 * @brief 	Here, we define the algorithm classes for erosion dynamics.
 * @details CK and SYCL version.
 * @author	Shuang Li and Xiangyu Hu
 */
#ifndef EROSION_DYNAMICS_2ND_H
#define EROSION_DYNAMICS_2ND_H

#include "continuum_integration_1st_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
// step2 inner
template <typename...>
class PlasticAcousticStepWithErosion2ndHalf;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStepWithErosion2ndHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStepWithErosion<Interaction<Inner<Parameters...>>>
{
    using PlasticKernel = typename PlasticContinuum::PlasticKernel;
    using BaseInteraction = PlasticAcousticStepWithErosion<Interaction<Inner<Parameters...>>>;

  public:
    explicit PlasticAcousticStepWithErosion2ndHalf(Inner<Parameters...> &inner_relation);
    virtual ~PlasticAcousticStepWithErosion2ndHalf() {};

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
        Real *rho_, *drho_dt_, *p_;
        Matd *velocity_gradient_;
        Mat3d *stress_tensor_3D_, *strain_tensor_3D_, *stress_rate_3D_, *strain_rate_3D_;
        Mat3d *total_stress_tensor_3D_, *viscous_stress_tensor_3D_, *shear_stress_tensor_3D_;
        Vecd *shear_vel_, *n_;
        Real *friction_angle_, *cohesion_, *reduction_para_, *test_;
        int *plastic_label_, *indicator_;
        Real particle_spacing_;
    };

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
    Real dv_particle_spacing_;
};

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStepWithErosion2ndHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStepWithErosion<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using BaseInteraction = PlasticAcousticStepWithErosion<Interaction<Contact<Parameters...>>>;

  public:
    explicit PlasticAcousticStepWithErosion2ndHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~PlasticAcousticStepWithErosion2ndHalf() {};

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

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStepWithErosion2ndHalf<Contact<Fluid, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStepWithErosion<Interaction<Contact<Parameters...>>>, public Interaction<Fluid>
{
    using BaseInteraction = PlasticAcousticStepWithErosion<Interaction<Contact<Parameters...>>>;

  public:
    explicit PlasticAcousticStepWithErosion2ndHalf(Contact<Parameters...> &fluid_contact_relation);
    virtual ~PlasticAcousticStepWithErosion2ndHalf() {};

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
        Vecd *shear_vel_; /*For erosion*/
        Real *fluid_Vol_, *fluid_p_;
        Vecd *fluid_vel_;
    };

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

using PlasticAcousticStepWithErosion2ndHalfWithWallRiemannCK =
    PlasticAcousticStepWithErosion2ndHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrection>,
                               Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrection>>;
} // namespace continuum_dynamics
} // namespace SPH
#endif // EROSION_DYNAMICS_2ND_H