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
 * @file 	repulsion_factor_summation.h
 * @brief 	TBD.
 * @details TBD.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef REPULSION_FACTOR_SUMMATION_H
#define REPULSION_FACTOR_SUMMATION_H

#include "base_contact_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{

template <typename... InteractionTypes>
class RepulsionFactorSummation;

template <class DataDelegationType>
class RepulsionFactorSummation<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    RepulsionFactorSummation(BaseRelationType &base_relation, const std::string &variable_name)
        : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
          repulsion_factor_(this->particles_->template registerStateVariable<Real>(variable_name)){};
    virtual ~RepulsionFactorSummation(){};

  protected:
    Real *repulsion_factor_;
};

template <>
class RepulsionFactorSummation<Inner<>> : public RepulsionFactorSummation<Base, DataDelegateInner>
{
  public:
    explicit RepulsionFactorSummation(SelfSurfaceContactRelation &self_contact_relation);
    virtual ~RepulsionFactorSummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real offset_W_ij_;
};
using SelfContactFactorSummation = RepulsionFactorSummation<Inner<>>;

template <>
class RepulsionFactorSummation<Contact<>> : public RepulsionFactorSummation<Base, DataDelegateContact>
{
  public:
    explicit RepulsionFactorSummation(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~RepulsionFactorSummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real> offset_W_ij_;
};
using ContactFactorSummation = RepulsionFactorSummation<Contact<>>;
/**
 * @class ShellContactFactor
 * @brief Computing the contact density due to shell contact using a
 * 		 surface integral being solved by Gauss-Legendre quadrature integration.
 *     This class can only be used when there's the source body only has contact to shell bodies,
 *     otherwise the contact density will be overwritten by that of solid contact bodies
 */
class ShellContactFactor : public RepulsionFactorSummation<Base, DataDelegateContact>
{
  public:
    explicit ShellContactFactor(ShellSurfaceContactRelation &solid_body_contact_relation);
    virtual ~ShellContactFactor(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    Kernel *kernel_;

    Real particle_spacing_;
    StdVec<Real> calibration_factor_;
    StdVec<Real> offset_W_ij_;
    StdVec<Real *> contact_Vol_;

    /** Abscissas and weights for Gauss-Legendre quadrature integration with n=3 nodes */
    const StdVec<Real> three_gaussian_points_ = {-0.7745966692414834, 0.0, 0.7745966692414834};
    const StdVec<Real> three_gaussian_weights_ = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
};

/**
 * @class ShellSelfContactFactorSummation
 * @brief Computing the contact density due to shell contact using dummy particle.
 */
class ShellSelfContactFactorSummation : public RepulsionFactorSummation<Base, DataDelegateInner>
{
  public:
    explicit ShellSelfContactFactorSummation(ShellSelfContactRelation &self_contact_relation);
    void interaction(size_t index_i, Real dt = 0.0);
};
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_FACTOR_SUMMATION_H
