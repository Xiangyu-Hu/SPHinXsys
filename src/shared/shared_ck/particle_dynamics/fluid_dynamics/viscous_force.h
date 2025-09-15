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
#include "viscosity.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... RelationTypes>
class ViscousForceCK;

template <typename ViscosityType, class KernelCorrectionType,
          template <typename...> class RelationType, typename... Parameters>
class ViscousForceCK<Base, ViscosityType, KernelCorrectionType, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
  public:
    template <class BaseRelationType>
    explicit ViscousForceCK(BaseRelationType &base_relation);
    virtual ~ViscousForceCK() {};

  protected:
    using ViscosityModel = ViscosityType;
    ViscosityType &viscosity_model_;
    KernelCorrectionType kernel_correction_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_viscous_force_;
    Real smoothing_length_sq_;
};

template <typename ViscosityType, class KernelCorrectionType, typename... Parameters>
class ViscousForceCK<Inner<WithUpdate, ViscosityType, KernelCorrectionType, Parameters...>>
    : public ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Inner<Parameters...>>,
      public ForcePriorCK
{
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using InterParticleViscosity = typename ViscosityType::InterParticleViscosity;
    using BaseViscousForceType = ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Inner<Parameters...>>;

  public:
    explicit ViscousForceCK(Inner<Parameters...> &inner_relation);
    virtual ~ViscousForceCK() {};

    class InteractKernel : public BaseViscousForceType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        InterParticleViscosity inter_particle_viscosity_;
        CorrectionKernel correction_;
        Real *Vol_;
        Vecd *vel_, *viscous_force_;
        Real smoothing_length_sq_;
    };
};

template <typename ViscosityType, class KernelCorrectionType, typename... Parameters>
class ViscousForceCK<Contact<Wall, ViscosityType, KernelCorrectionType, Parameters...>>
    : public ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Contact<Parameters...>>,
      public Interaction<Wall>
{
    using CorrectionKernel = typename KernelCorrectionType::ComputingKernel;
    using OneSideViscosity = typename ViscosityType::OneSideViscosity;
    using BaseViscousForceType = ViscousForceCK<Base, ViscosityType, KernelCorrectionType, Contact<Parameters...>>;

  public:
    explicit ViscousForceCK(Contact<Parameters...> &contact_relation);
    virtual ~ViscousForceCK() {};

    class InteractKernel : public BaseViscousForceType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        OneSideViscosity one_side_viscosity_;
        CorrectionKernel correction_;
        Real *Vol_, *contact_Vol_;
        Vecd *vel_, *viscous_force_, *wall_vel_ave_;
        Real smoothing_length_sq_;
    };

  protected:
    template <typename...>
    friend class FSI::ViscousForceFromFluid;
};

using ViscousForceInnerCK = ViscousForceCK<Inner<WithUpdate, Viscosity, NoKernelCorrectionCK>>;
using ViscousForceWithWallCK = ViscousForceCK<Inner<WithUpdate, Viscosity, NoKernelCorrectionCK>,
                                              Contact<Wall, Viscosity, NoKernelCorrectionCK>>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // VISCOUS_FORCE_H