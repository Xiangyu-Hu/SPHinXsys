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
 * @file 	erosion_dynamics_1st_ck.h
 * @brief 	Here, we define the algorithm classes for erosion dynamics.
 * @details CK and SYCL version.
 * @author	Shuang Li and Xiangyu Hu
 */
#ifndef EROSION_DYNAMICS_1ST_H
#define EROSION_DYNAMICS_1ST_H

#include "continuum_integration_1st_ck.h"

namespace SPH
{
namespace continuum_dynamics
{
template <class BaseInteractionType>
class PlasticAcousticStepWithErosion : public PlasticAcousticStep<BaseInteractionType>
{
public:
    template <class DynamicsIdentifier>
    explicit PlasticAcousticStepWithErosion(DynamicsIdentifier &identifier);
    virtual ~PlasticAcousticStepWithErosion() {};

protected:
    //Erosion related variables
    DiscreteVariable<Vecd> *dv_shear_vel_;          /* Fluid shear velocity */
    DiscreteVariable<Vecd> *dv_n_;          /* Normal direction */
    DiscreteVariable<Real> *dv_friction_angle_, *dv_cohesion_;     /* DP parameters */
    DiscreteVariable<Real> *dv_reduction_para_;      /* Strength reduction */
    DiscreteVariable<Real> *dv_test_;      /* Just for test */
    DiscreteVariable<int> *dv_plastic_label_;         /* Identify the eroded soil */
    DiscreteVariable<Mat3d> *dv_total_stress_tensor_3D_, *dv_viscous_stress_tensor_3D_, *dv_shear_stress_tensor_3D_;     
    DiscreteVariable<int> *dv_indicator_;
};


template <typename...>
class PlasticAcousticStepWithErosion1stHalf;

template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
class PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStepWithErosion<Interaction<Inner<Parameters...>>>
{
    using BaseInteraction = PlasticAcousticStepWithErosion<Interaction<Inner<Parameters...>>>;

  public:
    explicit PlasticAcousticStepWithErosion1stHalf(Inner<Parameters...> &inner_relation);
    virtual ~PlasticAcousticStepWithErosion1stHalf() {};

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
        Mat3d *stress_tensor_3D_, *viscous_stress_tensor_3D_, *shear_stress_tensor_3D_;
        int *plastic_label_;
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
class PlasticAcousticStepWithErosion1stHalf<Contact<Wall, RiemannSolverType, KernelCorrectionType, Parameters...>>
    : public PlasticAcousticStepWithErosion<Interaction<Contact<Parameters...>>>, public Interaction<Wall>
{
    using BaseInteraction = PlasticAcousticStepWithErosion<Interaction<Contact<Parameters...>>>;

  public:
    explicit PlasticAcousticStepWithErosion1stHalf(Contact<Parameters...> &wall_contact_relation);
    virtual ~PlasticAcousticStepWithErosion1stHalf() {};

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
        Mat3d *stress_tensor_3D_, *viscous_stress_tensor_3D_;

        Real *wall_Vol_;
        Vecd *wall_acc_ave_;
    };

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

using PlasticAcousticStepWithErosion1stHalfWithWallRiemannCK =
    PlasticAcousticStepWithErosion1stHalf<Inner<OneLevel, AcousticRiemannSolverCK, NoKernelCorrection>,
                               Contact<Wall, AcousticRiemannSolverCK, NoKernelCorrection>>;
} // namespace continuum_dynamics
} // namespace SPH
#endif // EROSION_DYNAMICS_1ST_H