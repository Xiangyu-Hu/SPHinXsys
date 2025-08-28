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
 * @file 	contact_repulsion.h
 * @brief 	TBD.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONTACT_REPULSION_H
#define CONTACT_REPULSION_H

#include "base_contact_dynamics.h"
#include "force_prior.hpp"

namespace SPH
{
namespace solid_dynamics
{

template <typename... InteractionTypes>
class RepulsionForce;

template <class DataDelegationType>
class RepulsionForce<Base, DataDelegationType>
    : public ForcePrior, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    RepulsionForce(BaseRelationType &base_relation, const std::string &variable_name)
        : ForcePrior(base_relation.getSPHBody(), variable_name), DataDelegationType(base_relation),
          repulsion_force_(this->particles_->template registerStateVariable<Vecd>(variable_name)),
          Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")){};
    virtual ~RepulsionForce(){};

  protected:
    Vecd *repulsion_force_;
    Real *Vol_;
};

template <>
class RepulsionForce<Contact<Inner<>>> : public RepulsionForce<Base, DataDelegateInner>
{
  public:
    explicit RepulsionForce(BaseInnerRelation &self_contact_relation);
    virtual ~RepulsionForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    Real *self_repulsion_factor_;
    Vecd *vel_;
    Real contact_impedance_;
};
using SelfContactForce = RepulsionForce<Contact<Inner<>>>;

template <>
class RepulsionForce<Contact<>> : public RepulsionForce<Base, DataDelegateContact>
{
  public:
    explicit RepulsionForce(BaseContactRelation &solid_body_contact_relation);
    virtual ~RepulsionForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    Real *repulsion_factor_;
    StdVec<Solid *> contact_solids_;
    StdVec<Real> contact_stiffness_ave_;
    StdVec<Real *> contact_repulsion_factor_, contact_Vol_;
};
using ContactForce = RepulsionForce<Contact<>>;

template <> // Computing the repulsion force from a rigid wall.
class RepulsionForce<Contact<Wall>> : public RepulsionForce<Base, DataDelegateContact>
{
  public:
    explicit RepulsionForce(BaseContactRelation &solid_body_contact_relation);
    virtual ~RepulsionForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    Real *repulsion_factor_;
    StdVec<Real *> contact_Vol_;
};
using ContactForceFromWall = RepulsionForce<Contact<Wall>>;

template <> // Computing the repulsion force acting on a rigid wall.
class RepulsionForce<Wall, Contact<>> : public RepulsionForce<Base, DataDelegateContact>
{
  public:
    explicit RepulsionForce(BaseContactRelation &solid_body_contact_relation);
    virtual ~RepulsionForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Solid *> contact_solids_;
    StdVec<Real *> contact_repulsion_factor_, contact_Vol_;
};
using ContactForceToWall = RepulsionForce<Wall, Contact<>>;
} // namespace solid_dynamics
} // namespace SPH
#endif // CONTACT_REPULSION_H
