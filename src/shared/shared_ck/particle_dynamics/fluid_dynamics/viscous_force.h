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

namespace SPH
{
namespace fluid_dynamics
{
template <typename...>
class Viscosity;

template <>
class Viscosity<Constant>
{
  public:
    Viscosity(BaseParticles *particles1, BaseParticles *particles2)
        : mu1_(DynamicCast<Fluid>(this, particles1->getBaseMaterial()).ReferenceViscosity()),
          mu2_(DynamicCast<Fluid>(this, particles2->getBaseMaterial()).ReferenceViscosity()){};

    class ComputingKernel : public PairGeomAverageFixed<Real>
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy, Viscosity<Constant> &encloser,
                        ComputingKernelType &computing_kernel)
            : PairGeomAverageFixed<Real>(encloser.mu1_, encloser.mu2_){};
    };

  protected:
    Real mu1_, mu2_;
};

template <>
class Viscosity<Variable>
{
  public:
    Viscosity(BaseParticles *particles1, BaseParticles *particles2)
        : dv_mu1_(particles1->getVariableByName<Real>("VariableViscosity")),
          dv_mu2_(particles2->getVariableByName<Real>("VariableViscosity")){};

    class ComputingKernel : public PairGeomAverageVariable<Real>
    {
      public:
        template <class ExecutionPolicy, class ComputingKernelType>
        ComputingKernel(const ExecutionPolicy &ex_policy, Viscosity<Variable> &encloser,
                        ComputingKernelType &computing_kernel)
            : PairGeomAverageVariable<Real>(
                  encloser.dv_mu1_->DelegatedDataField(ex_policy),
                  encloser.dv_mu2_->DelegatedDataField(ex_policy)){};
    };

  protected:
    DiscreteVariable<Real> *dv_mu1_, *dv_mu2_;
};

template <typename... RelationTypes>
class ViscousForceCK;

template <template <typename...> class RelationType, typename... Parameters>
class ViscousForceCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>, public ForcePriorCK
{
  public:
    template <class BaseRelationType>
    explicit ViscousForceCK(BaseRelationType &base_relation);
    virtual ~ViscousForceCK(){};

    class InteractKernel
        : public Interaction<RelationType<Parameters...>>::InteractKernel,
          public ForcePriorCK::UpdateKernel
    {
      public:
        template <class ExecutionPolicy, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy,
                       ViscousForceCK<Base, RelationType<Parameters...>> &encloser,
                       Args &&...args);

      protected:
        Real *rho_, *mass_, *Vol_;
        Vecd *vel_, *viscous_force_;
        Real smoothing_length_;
    };

  protected:
    DiscreteVariable<Real> *dv_rho_, *dv_mass_, *dv_Vol_;
    DiscreteVariable<Vecd> *dv_vel_, *dv_viscous_force_;
    Real smoothing_length_;
};

} // namespace fluid_dynamics
} // namespace SPH
#endif // VISCOUS_FORCE_H