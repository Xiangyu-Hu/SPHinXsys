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
 * @file 	repulsion_force.h
 * @brief 	TBD.
 * @details TBD.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef REPULSION_FORCE_H
#define REPULSION_FORCE_H

#include "force_prior_ck.h"

namespace SPH
{
namespace solid_dynamics
{

template <typename...>
class RepulsionForceCK;

template <typename... Parameters>
class RepulsionForceCK<Base, Contact<Parameters...>>
    : public Interaction<Contact<Parameters...>>, public ForcePriorCK
{
  public:
    explicit RepulsionForceCK(Contact<Parameters...> &contact_relation, Real numerical_damping = 0.0);
    virtual ~RepulsionForceCK() {};

  protected:
    SolidContact &solid_contact_;
    Real numerical_damping_;
    Real stiffness_, impedance_;
    DiscreteVariable<Real> *dv_repulsion_factor_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_vel_;
    DiscreteVariable<Vecd> *dv_repulsion_force_;
};

template <typename... Parameters>
class RepulsionForceCK<Contact<WithUpdate, Parameters...>>
    : public RepulsionForceCK<Base, Contact<Parameters...>>
{
    using BaseInteractionType = RepulsionForceCK<Base, Contact<Parameters...>>;

  public:
    template <typename... Args>
    RepulsionForceCK(Contact<Parameters...> &contact_relation, Args &&...args);
    virtual ~RepulsionForceCK() {};

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real numerical_damping_;
        Real stiffness_ave_, contact_stiffness_, impedance_, contact_impedance_;
        Real *repulsion_factor_, *contact_repulsion_factor_;
        Real *Vol_, *contact_Vol_;
        Vecd *vel_, *contact_vel_;
        Vecd *n_, *contact_n_;
        Vecd *repulsion_force_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_n_;
    StdVec<Real> contact_stiffness_, contact_impedance_;
    StdVec<DiscreteVariable<Real> *> dv_contact_repulsion_factor_;
    StdVec<DiscreteVariable<Vecd> *> dv_contact_vel_, dv_contact_n_;
};

template <typename... Parameters>
class RepulsionForceCK<Contact<WithUpdate, Wall, Parameters...>>
    : public RepulsionForceCK<Base, Contact<Parameters...>>, public Interaction<Wall>
{
    using BaseInteractionType = RepulsionForceCK<Base, Contact<Parameters...>>;

  public:
    template <typename... Args>
    RepulsionForceCK(Contact<Parameters...> &contact_relation, Args &&...args);
    virtual ~RepulsionForceCK() {};

    class InteractKernel : public BaseInteractionType::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);
        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real numerical_damping_;
        Real stiffness_, impedance_;
        Real *repulsion_factor_;
        Real *Vol_, *contact_Vol_;
        Vecd *vel_, *contact_vel_;
        Vecd *contact_n_;
        Vecd *repulsion_force_;
    };
};
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_FORCE_H