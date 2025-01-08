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
 * @file 	viscous_force.h
 * @brief Here, we define the algorithm classes for computing viscous forces in fluids.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef VISCOUS_FORCE_H
#define VISCOUS_FORCE_H

#include "base_fluid_dynamics.h"
#include "force_prior_ck.h"
#include "interaction_ck.hpp"
#include "kernel_correction_ck.hpp"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... RelationTypes>
class ViscousForceCK;

template <typename ViscosityType, class KernelCorrectionType,
          template <typename...> class RelationType, typename... Parameters>
class ViscousForceCK<Base, ViscosityType, KernelCorrectionType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>, public ForcePriorCK
{
    using ViscosityKernel = typename ViscosityType::ComputingKernel;
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;

  public:
    template <class BaseRelationType>
    explicit ViscousForceCK(BaseRelationType &base_relation);
    virtual ~ViscousForceCK() {};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args);

      protected:
        ViscosityKernel viscosity_;
        CorrectionKernel correction_;
        Real *Vol_;
        Vecd *vel_, *viscous_force_;
        Real smoothing_length_sq_;
    };

  protected:
    template <typename...>
    friend class FSI::ViscousForceFromFluid;

    ViscosityType viscosity_method_;
    KernelCorrectionType kernel_correction_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_viscous_force_;
    Real smoothing_length_sq_;
};

template <typename ViscosityType, class KernelCorrectionType, typename... Parameters>
class ViscousForceCK<Inner<WithUpdate, ViscosityType, KernelCorrectionType, Parameters...>>
    : public ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Inner<Parameters...>>
{
    using BaseViscousForceType = ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Inner<Parameters...>>;

  public:
    explicit ViscousForceCK(Relation<Inner<Parameters...>> &inner_relation)
        : BaseViscousForceType(inner_relation) {};
    virtual ~ViscousForceCK() {};

    class InteractKernel : public BaseViscousForceType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseViscousForceType::InteractKernel(ex_policy, encloser){};
        void interact(size_t index_i, Real dt = 0.0);
    };
};

template <typename ViscosityType, class KernelCorrectionType, typename... Parameters>
class ViscousForceCK<Contact<Wall, ViscosityType, KernelCorrectionType, Parameters...>>
    : public ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Contact<Wall, Parameters...>>
{
    using BaseViscousForceType = ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Contact<Wall, Parameters...>>;

  public:
    explicit ViscousForceCK(Relation<Contact<Parameters...>> &contact_relation)
        : BaseViscousForceType(contact_relation) {};
    virtual ~ViscousForceCK() {};

    class InteractKernel : public BaseViscousForceType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index)
            : BaseViscousForceType::InteractKernel(ex_policy, encloser, contact_index),
              wall_Vol_(encloser.dv_wall_Vol_[contact_index]->DelegatedData(ex_policy)),
              wall_vel_ave_(encloser.dv_wall_vel_ave_[contact_index]->DelegatedData(ex_policy)){};
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *wall_Vol_;
        Vecd *wall_vel_ave_;
    };
};

template <typename...>
class Viscosity;

template <>
class Viscosity<Constant>
{
  public:
    Viscosity(BaseParticles *particles)
        : mu_(DynamicCast<Fluid>(this, particles->getBaseMaterial()).ReferenceViscosity()) {};

    class ComputingKernel : public ParameterFixed<Real>
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : ParameterFixed<Real>(encloser.mu_){};
    };

  protected:
    Real mu_;
};

template <>
class Viscosity<Variable>
{
  public:
    Viscosity(BaseParticles *particles)
        : dv_mu_(particles->getVariableByName<Real>("VariableViscosity")) {};

    class ComputingKernel : public ParameterVariable<Real>
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : ParameterVariable<Real>(encloser.dv_mu_->DelegatedData(ex_policy)){};
    };

  protected:
    DiscreteVariable<Real> *dv_mu_;
};

using ViscousForceInnerCK = ViscousForceCK<Inner<WithUpdate, Viscosity<Constant>, NoKernelCorrectionCK>>;
using ViscousForceWithWallCK = ViscousForceCK<Inner<WithUpdate, Viscosity<Constant>, NoKernelCorrectionCK>,
                                              Contact<Wall, Viscosity<Constant>, NoKernelCorrectionCK>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // VISCOUS_FORCE_H