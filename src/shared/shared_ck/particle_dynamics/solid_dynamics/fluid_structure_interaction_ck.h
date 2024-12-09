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
 * @file 	fluid_structure_interaction_ck.h
 * @brief 	Here, we define the algorithm classes for fluid structure interaction.
 * @author	Xiangyu Hu
 */

#ifndef FLUID_STRUCTURE_INTERACTION_CK_H
#define FLUID_STRUCTURE_INTERACTION_CK_H

#include "base_material.h"
#include "force_prior_ck.hpp"
#include "interaction_ck.hpp"
#include "riemann_solver.h"

namespace SPH
{
namespace FSI
{
template <class KernelCorrectionType, typename... Parameters>
class ForceFromFluidCK : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    explicit ForceFromFluidCK(BaseContactRelation &contact_relation, const std::string &force_name);
    virtual ~ForceFromFluidCK(){};

    class InteractKernel
        : public Interaction<Contact<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);

      protected:
        Real *Vol_;
        Vecd *force_from_fluid_, *vel_ave_;
        CorrectionKernel contact_correction_;
        Real *contact_Vol_;
        Vecd *contact_vel_;
    };

  protected:
    Solid &solid_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_force_from_fluid_, *dv_vel_ave_;

    StdVec<Fluid *> contact_fluid_;
    StdVec<KernelCorrectionType> contact_kernel_correction_;
    StdVec<DiscreteVariable<Real> *> contact_Vol_;
    StdVec<DiscreteVariable<Vecd> *> contact_vel_;
};

template <typename...>
class ViscousForceFromFluidCK;

template <typename ViscousForceType, typename... Parameters>
class ViscousForceFromFluidCK<Contact<WithUpdate, ViscousForceType, Parameters...>>
    : public ForceFromFluidCK<decltype(ViscousForceType::kernel_correction_), Parameters...>
{

    using ViscosityType = decltype(ViscousForceType::viscosity_method_);
    using ViscosityKernel = typename ViscosityType::ComputingKernel;
    using BaseForceFromFluid = ForceFromFluidCK<decltype(ViscousForceType::kernel_correction_), Parameters...>;

  public:
    explicit ViscousForceFromFluidCK(BaseContactRelation &contact_relation);
    virtual ~ViscousForceFromFluidCK(){};
    class InteractKernel : public BaseForceFromFluid::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        ViscosityKernel viscosity_;
        Real smoothing_length_sq_;
    };

  protected:
    StdVec<ViscosityType> contact_viscosity_method_;
    StdVec<Real> contact_smoothing_length_sq_;
};
template <typename ViscousForceType>
using ViscousForceOnStructure = ViscousForceFromFluidCK<Contact<WithUpdate, ViscousForceType>>;

template <typename...>
class PressureForceFromFluidCK;

template <class AcousticStep2ndHalfType, typename... Parameters>
class PressureForceFromFluidCK<Contact<WithUpdate, AcousticStep2ndHalfType, Parameters...>>
    : public ForceFromFluidCK<decltype(AcousticStep2ndHalfType::kernel_correction_), Parameters...>
{
    using RiemannSolverType = decltype(AcousticStep2ndHalfType::riemann_solver_);
    using KernelCorrectionType = decltype(AcousticStep2ndHalfType::kernel_correction_);
    using BaseForceFromFluid = ForceFromFluidCK<decltype(AcousticStep2ndHalfType::kernel_correction_), Parameters...>;

  public:
    explicit PressureForceFromFluidCK(BaseContactRelation &contact_relation);
    virtual ~PressureForceFromFluidCK(){};

    class InteractKernel : public BaseForceFromFluid::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, UnsignedInt contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *acc_ave_, *n_;
        RiemannSolverType riemann_solver_;
        Real *contact_rho_, *contact_mass_, *contact_p_;
        Vecd *contact_force_prior_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_acc_ave_, *dv_n_;
    StdVec<RiemannSolverType> contact_riemann_solver_;
    StdVec<DiscreteVariable<Real> *> dv_contact_rho_, dv_contact_mass_, dv_contact_p_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_force_prior_;
};
template <typename AcousticStep2ndHalfType>
using PressureForceOnStructure = PressureForceFromFluidCK<Contact<WithUpdate, AcousticStep2ndHalfType>>;
} // namespace FSI
} // namespace SPH
#endif // FLUID_STRUCTURE_INTERACTION_CK_H