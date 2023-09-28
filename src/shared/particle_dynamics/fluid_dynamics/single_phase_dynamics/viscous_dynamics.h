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
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef VISCOUS_DYNAMICS_H
#define VISCOUS_DYNAMICS_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class BaseViscousAcceleration
 * @brief Base class for the viscosity force induced acceleration
 */
template <class DataDelegationType>
class BaseViscousAcceleration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseViscousAcceleration(BaseRelationType &base_relation);
    virtual ~BaseViscousAcceleration(){};

  protected:
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &vel_, &acc_prior_;
    Real mu_;
    Real smoothing_length_;
};

/**
 * @class ViscousAccelerationInner
 * @brief  the viscosity force induced acceleration
 */
class ViscousAccelerationInner : public BaseViscousAcceleration<FluidDataInner>
{
  public:
    explicit ViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAcceleration<FluidDataInner>(inner_relation){};
    virtual ~ViscousAccelerationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class AngularConservativeViscousAccelerationInner
 * @brief the viscosity force induced acceleration, a formulation for conserving
 * angular momentum, to be tested for its practical applications.
 */
class AngularConservativeViscousAccelerationInner : public BaseViscousAcceleration<FluidDataInner>
{
  public:
    explicit AngularConservativeViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAcceleration<FluidDataInner>(inner_relation){};
    virtual ~AngularConservativeViscousAccelerationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ViscousWallBoundary
 * @brief TBD
 */
class ViscousWallBoundary : public InteractionWithWall<BaseViscousAcceleration>
{
  public:
    ViscousWallBoundary(BaseContactRelation &wall_contact_relation)
        : InteractionWithWall<BaseViscousAcceleration>(wall_contact_relation){};
    virtual ~ViscousWallBoundary(){};
    void interaction(size_t index_i, Real dt = 0.0);
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

class ViscousAccelerationWithWall
    : public ComplexInteraction<ViscousAccelerationInner, ViscousWallBoundary>
{
  public:
    explicit ViscousAccelerationWithWall(ComplexRelation &fluid_wall_relation)
        : ComplexInteraction<ViscousAccelerationInner, ViscousWallBoundary>(
              fluid_wall_relation.getInnerRelation(), fluid_wall_relation.getContactRelation()){};
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // VISCOUS_DYNAMICS_H