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
 * @file 	force_on_structure.h
 * @brief 	Here, we define the algorithm classes for fluid structure interaction.
 * @author	Xiangyu Hu
 */

#ifndef FORCE_ON_STRUCTURE_H
#define FORCE_ON_STRUCTURE_H

#include "base_material.h"
#include "force_prior_ck.hpp"
#include "interaction_ck.hpp"
#include "riemann_solver_ck.hpp"
#include "viscosity.h"

namespace SPH
{
namespace FSI
{
template <class KernelCorrectionType, typename... Parameters>
class ForceFromFluid : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    template <class ContactRelationType>
    explicit ForceFromFluid(ContactRelationType &contact_relation, const std::string &force_name);
    virtual ~ForceFromFluid(){};

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

    StdVec<KernelCorrectionType> contact_kernel_correction_;
    StdVec<DiscreteVariable<Real> *> contact_Vol_;
    StdVec<DiscreteVariable<Vecd> *> contact_vel_;
};

template <typename...>
class ViscousForceFromFluid;

template <typename ViscousForceType, typename... Parameters>
class ViscousForceFromFluid<Contact<WithUpdate, ViscousForceType, Parameters...>>
    : public ForceFromFluid<decltype(ViscousForceType::kernel_correction_), Parameters...>
{

    using ViscosityType = typename ViscousForceType::ViscosityModel;
    using ViscosityKernel = typename ViscosityType::ComputingKernel;
    using BaseForceFromFluid = ForceFromFluid<decltype(ViscousForceType::kernel_correction_), Parameters...>;

  public:
    template <class ContactRelationType>
    explicit ViscousForceFromFluid(ContactRelationType &contact_relation);
    virtual ~ViscousForceFromFluid(){};
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
    StdVec<ViscosityType *> contact_viscosity_model_;
    StdVec<Real> contact_smoothing_length_sq_;
};
template <typename ViscousForceType>
using ViscousForceOnStructure = ViscousForceFromFluid<Contact<WithUpdate, ViscousForceType>>;

template <typename...>
class PressureForceFromFluid;

template <class AcousticStep2ndHalfType, typename... Parameters>
class PressureForceFromFluid<Contact<WithUpdate, AcousticStep2ndHalfType, Parameters...>>
    : public ForceFromFluid<decltype(AcousticStep2ndHalfType::kernel_correction_), Parameters...>
{
    using RiemannSolverType = decltype(AcousticStep2ndHalfType::riemann_solver_);
    using FluidType = typename RiemannSolverType::SourceFluid;
    using KernelCorrectionType = decltype(AcousticStep2ndHalfType::kernel_correction_);
    using BaseForceFromFluid = ForceFromFluid<decltype(AcousticStep2ndHalfType::kernel_correction_), Parameters...>;

  public:
    template <class ContactRelationType>
    explicit PressureForceFromFluid(ContactRelationType &contact_relation);
    virtual ~PressureForceFromFluid(){};

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
using PressureForceOnStructure = PressureForceFromFluid<Contact<WithUpdate, AcousticStep2ndHalfType>>;
} // namespace FSI
} // namespace SPH
#endif // FORCE_ON_STRUCTURE_H