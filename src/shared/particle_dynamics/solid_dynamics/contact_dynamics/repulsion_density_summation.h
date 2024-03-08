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
 * @file 	repulsion_density_summation.h
 * @brief 	TBD.
 * @details TBD.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef REPULSION_DENSITY_SUMMATION_H
#define REPULSION_DENSITY_SUMMATION_H

#include "base_contact_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
typedef DataDelegateContact<SolidParticles, SolidParticles> ContactDynamicsData;
typedef DataDelegateContact<SolidParticles, SolidParticles> ContactWithWallData;

template <typename... InteractionTypes>
class RepulsionDensitySummation;

template <class DataDelegationType>
class RepulsionDensitySummation<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    RepulsionDensitySummation(BaseRelationType &base_relation, const std::string &variable_name)
        : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
          repulsion_density_(*this->particles_->template registerSharedVariable<Real>(variable_name)){};
    virtual ~RepulsionDensitySummation(){};

  protected:
    StdLargeVec<Real> &repulsion_density_;
};

template <>
class RepulsionDensitySummation<Inner<>> : public RepulsionDensitySummation<Base, SolidDataInner>
{
  public:
    explicit RepulsionDensitySummation(BaseInnerRelation &self_contact_relation);
    virtual ~RepulsionDensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &mass_;
    Real offset_W_ij_;
};
using SelfContactDensitySummation = RepulsionDensitySummation<Inner<>>;

template <>
class RepulsionDensitySummation<Contact<>> : public RepulsionDensitySummation<Base, ContactDynamicsData>
{
  public:
    explicit RepulsionDensitySummation(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~RepulsionDensitySummation(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    UniquePtrKeeper<Kernel> kernel_keeper_;
    StdLargeVec<Real> &mass_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    StdVec<Real> offset_W_ij_;
};
using ContactDensitySummation = RepulsionDensitySummation<Contact<>>;
/**
 * @class ShellContactDensity
 * @brief Computing the contact density due to shell contact using a
 * 		 surface integral being solved by Gauss-Legendre quadrature integration.
 */
class ShellContactDensity : public RepulsionDensitySummation<Base, ContactDynamicsData>
{
  public:
    explicit ShellContactDensity(SurfaceContactRelation &solid_body_contact_relation);
    virtual ~ShellContactDensity(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Solid &solid_;
    Kernel *kernel_;

    Real particle_spacing_;
    StdVec<Real> calibration_factor_;
    StdVec<Real> offset_W_ij_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;

    /** Abscissas and weights for Gauss-Legendre quadrature integration with n=3 nodes */
    const StdVec<Real> three_gaussian_points_ = {-0.7745966692414834, 0.0, 0.7745966692414834};
    const StdVec<Real> three_gaussian_weights_ = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
};

/**
 * @class ShellSelfContactDensityUsingDummyParticles
 * @brief Computing the contact density due to shell contact using dummy particle.
 */
class ShellSelfContactDensityUsingDummyParticles : public RepulsionDensitySummation<Inner<>>
{
  public:
    explicit ShellSelfContactDensityUsingDummyParticles(ShellSelfContactRelation &self_contact_relation);

  private:
    UniquePtrKeeper<Kernel> kernel_keeper_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif // REPULSION_DENSITY_SUMMATION_H
