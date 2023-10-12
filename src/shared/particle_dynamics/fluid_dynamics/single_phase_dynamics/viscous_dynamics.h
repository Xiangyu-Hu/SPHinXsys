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

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class ViscousAcceleration;

template <class DataDelegationType>
class ViscousAcceleration<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit ViscousAcceleration(BaseRelationType &base_relation);
    virtual ~ViscousAcceleration(){};

  protected:
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &vel_, &acc_prior_;
    Real mu_;
    Real smoothing_length_;
};

template <>
class ViscousAcceleration<Inner<>> : public ViscousAcceleration<FluidDataInner>
{
  public:
    explicit ViscousAcceleration(BaseInnerRelation &inner_relation)
        : ViscousAcceleration<FluidDataInner>(inner_relation){};
    virtual ~ViscousAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using ViscousAccelerationInner = ViscousAcceleration<Inner<>>;

template <>
class ViscousAcceleration<AngularConservative<Inner<>>>
    : public ViscousAcceleration<FluidDataInner>
{
  public:
    explicit ViscousAcceleration(BaseInnerRelation &inner_relation)
        : ViscousAcceleration<FluidDataInner>(inner_relation){};
    virtual ~ViscousAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

template <>
class ViscousAcceleration<ContactWall<>> : public InteractionWithWall<ViscousAcceleration>
{
  public:
    explicit ViscousAcceleration(BaseContactRelation &wall_contact_relation)
        : InteractionWithWall<ViscousAcceleration>(wall_contact_relation){};
    virtual ~ViscousAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

class ViscousAccelerationWithWall
    : public ComplexInteraction<ViscousAcceleration<Inner<>, ContactWall<>>>
{
  public:
    explicit ViscousAccelerationWithWall(ComplexRelation &fluid_wall_relation)
        : ComplexInteraction<ViscousAcceleration<Inner<>, ContactWall<>>>(
              fluid_wall_relation.getInnerRelation(), fluid_wall_relation.getContactRelation()){};
};

template <>
class ViscousAcceleration<Contact<>> : public ViscousAcceleration<FluidContactData>
{
  public:
    explicit ViscousAcceleration(BaseContactRelation &contact_relation);
    virtual ~ViscousAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real> contact_mu_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};

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