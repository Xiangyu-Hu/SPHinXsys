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
 * @file 	viscous_dynamics.h
 * @brief Here, we define the algorithm classes for computing viscous forces in fluids.
 * @details TBD.
 * @author	Xiangyu Hu
 */

#ifndef VISCOUS_DYNAMICS_H
#define VISCOUS_DYNAMICS_H

#include "base_fluid_dynamics.h"
#include "force_prior.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class ViscousForce;

template <class DataDelegationType>
class ViscousForce<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit ViscousForce(BaseRelationType &base_relation);
    virtual ~ViscousForce(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_;
    StdLargeVec<Vecd> &vel_, &viscous_force_;
    Real mu_;
    Real smoothing_length_;
};

template <>
class ViscousAcceleration<Inner<>>
    : public ViscousAcceleration<FluidDataInner>, public ForcePrior
{
  public:
    explicit ViscousAcceleration(BaseInnerRelation &inner_relation)
        : ViscousAcceleration<FluidDataInner>(inner_relation),
          ForcePrior(&base_particles_, "ViscousForce"){};
    virtual ~ViscousAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using ViscousForceInner = ViscousForce<Inner<>>;

template <>
class ViscousAcceleration<AngularConservative<Inner<>>>
    : public ViscousAcceleration<FluidDataInner>, public ForcePrior
{
  public:
    explicit ViscousAcceleration(BaseInnerRelation &inner_relation)
        : ViscousAcceleration<FluidDataInner>(inner_relation),
          ForcePrior(&base_particles_, "ViscousForce"){};
    virtual ~ViscousAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

using BaseViscousForceWithWall = InteractionWithWall<ViscousForce>;
template <>
class ViscousForce<Contact<Wall>> : public BaseViscousForceWithWall
{
  public:
    explicit ViscousForce(BaseContactRelation &wall_contact_relation)
        : BaseViscousForceWithWall(wall_contact_relation){};
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <>
class ViscousForce<Contact<>> : public ViscousForce<FluidContactData>
{
  public:
    explicit ViscousForce(BaseContactRelation &contact_relation);
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real> contact_mu_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};

using ViscousForceWithWall = ComplexInteraction<ViscousForce<Inner<>, Contact<Wall>>>;

/**
 * @class VorticityInner
 * @brief  compute vorticity in the fluid field
 */
class VorticityInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit VorticityInner(BaseInnerRelation &inner_relation);
    virtual ~VorticityInner(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<AngularVecd> vorticity_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // VISCOUS_DYNAMICS_H