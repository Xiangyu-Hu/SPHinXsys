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

class FixedViscosity : public PairGeomAverageFixed<Real>
{
  public:
    FixedViscosity(BaseParticles *particles1, BaseParticles *particles2)
        : PairGeomAverageFixed<Real>(
              DynamicCast<Fluid>(this, particles1->getBaseMaterial()).ReferenceViscosity(),
              DynamicCast<Fluid>(this, particles2->getBaseMaterial()).ReferenceViscosity()){};
    explicit FixedViscosity(BaseParticles *particles)
        : FixedViscosity(particles, particles){};
    virtual ~FixedViscosity(){};
};

class VariableViscosity : public PairGeomAverageVariable<Real>
{
  public:
    VariableViscosity(BaseParticles *particles1, BaseParticles *particles2)
        : PairGeomAverageVariable<Real>(
              *particles1->getVariableDataByName<Real>("VariableViscosity"),
              *particles2->getVariableDataByName<Real>("VariableViscosity")){};
    explicit VariableViscosity(BaseParticles *particles)
        : VariableViscosity(particles, particles){};
    virtual ~VariableViscosity(){};
};

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
    StdLargeVec<Real> &rho_, &mass_, &Vol_;
    StdLargeVec<Vecd> &vel_, &viscous_force_;
    Real smoothing_length_;
};

template <typename ViscosityType>
class ViscousForce<Inner<>, ViscosityType>
    : public ViscousForce<DataDelegateInner>, public ForcePrior
{
  public:
    explicit ViscousForce(BaseInnerRelation &inner_relation);
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    ViscosityType mu_;
};
using ViscousForceInner = ViscousForce<Inner<>, FixedViscosity>;

template <typename ViscosityType>
class ViscousForce<Inner<AngularConservative>, ViscosityType>
    : public ViscousForce<DataDelegateInner>, public ForcePrior
{
  public:
    explicit ViscousForce(BaseInnerRelation &inner_relation)
        : ViscousForce<DataDelegateInner>(inner_relation),
          ForcePrior(particles_, "ViscousForce"),
          mu_(particles_){};
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    ViscosityType mu_;
};

using BaseViscousForceWithWall = InteractionWithWall<ViscousForce>;

template <typename ViscosityType>
class ViscousForce<Contact<Wall>, ViscosityType> : public BaseViscousForceWithWall
{
  public:
    explicit ViscousForce(BaseContactRelation &wall_contact_relation)
        : BaseViscousForceWithWall(wall_contact_relation),
          mu_(particles_){};
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    ViscosityType mu_;
};

template <typename ViscosityType>
class ViscousForce<Contact<Wall, AngularConservative>, ViscosityType> : public BaseViscousForceWithWall
{
  public:
    explicit ViscousForce(BaseContactRelation &wall_contact_relation);
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    ViscosityType mu_;
    StdLargeVec<Vecd> &distance_from_wall_;
};

template <typename ViscosityType>
class ViscousForce<Contact<>, ViscosityType> : public ViscousForce<DataDelegateContact>
{
  public:
    explicit ViscousForce(BaseContactRelation &contact_relation);
    virtual ~ViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<ViscosityType> contact_mu_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
    StdVec<StdLargeVec<Real> *> wall_Vol_;
};

using ViscousForceWithWall = ComplexInteraction<ViscousForce<Inner<>, Contact<Wall>>, FixedViscosity>;
using MultiPhaseViscousForceWithWall = ComplexInteraction<ViscousForce<Inner<>, Contact<>, Contact<Wall>>, FixedViscosity>;

template <typename... FormulationType>
using NonNewtonianViscousForceWithWall =
    ComplexInteraction<ViscousForce<Inner<FormulationType...>, Contact<Wall, FormulationType...>>, VariableViscosity>;

/**
 * @class VorticityInner
 * @brief  compute vorticity in the fluid field
 */
class VorticityInner : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VorticityInner(BaseInnerRelation &inner_relation);
    virtual ~VorticityInner(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &Vol_;
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<AngularVecd> &vorticity_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // VISCOUS_DYNAMICS_H