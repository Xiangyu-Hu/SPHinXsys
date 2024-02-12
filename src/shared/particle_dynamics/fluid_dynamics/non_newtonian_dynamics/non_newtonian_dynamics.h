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
 * @file non_newtonian_dynamics.h
 * @brief Here, we define the time integration algorithm classes for non-newtonian fluids.
 * @details We consider here weakly compressible fluids.
 * @author	Xiangyu Hu
 */

#ifndef NON_NEWTONIAN_DYNAMICS_H
#define NON_NEWTONIAN_DYNAMICS_H

#include "base_fluid_dynamics.h"
#include "fluid_integration.hpp"
#include "force_prior.h"

namespace SPH
{
namespace fluid_dynamics
{
template <typename... InteractionTypes>
class Oldroyd_BIntegration1stHalf;

template <>
class Oldroyd_BIntegration1stHalf<Inner<>> : public Integration1stHalfInnerRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> tau_, dtau_dt_;
};

using Integration1stHalfContactWallRiemann =
    Integration1stHalf<Contact<Wall>, AcousticRiemannSolver, NoKernelCorrection>;

template <>
class Oldroyd_BIntegration1stHalf<Contact<Wall>> : public Integration1stHalfContactWallRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &tau_;
};

template <typename... InteractionTypes>
class Oldroyd_BIntegration2ndHalf;

template <>
class Oldroyd_BIntegration2ndHalf<Inner<>> : public Integration2ndHalfInnerRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};

using Integration2ndHalfContactWallRiemann =
    Integration2ndHalf<Contact<Wall>, AcousticRiemannSolver>;

template <>
class Oldroyd_BIntegration2ndHalf<Contact<Wall>> : public Integration2ndHalfContactWallRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};

using Oldroyd_BIntegration1stHalfWithWall = ComplexInteraction<Oldroyd_BIntegration1stHalf<Inner<>, Contact<Wall>>>;
using Oldroyd_BIntegration2ndHalfWithWall = ComplexInteraction<Oldroyd_BIntegration2ndHalf<Inner<>, Contact<Wall>>>;

/**
 * @class GenearlizedNewtonianViscousForce
 * @brief Calculates the viscous force based on a generalized Newtonian fluid model
 * @note This needs the VelocityGradient & ShearRateDependentViscosity to work
 */
template <typename... InteractionTypes>
class GeneralizedNewtonianViscousForce;

template <class DataDelegationType>
class GeneralizedNewtonianViscousForce<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit GeneralizedNewtonianViscousForce(BaseRelationType &base_relation);
    virtual ~GeneralizedNewtonianViscousForce(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_;
    StdLargeVec<Vecd> &vel_, &viscous_force_;
    GeneralizedNewtonianFluid &generalized_newtonian_fluid_;
    Real smoothing_length_;
};

template <>
class GeneralizedNewtonianViscousForce<Inner<>> : public GeneralizedNewtonianViscousForce<FluidDataInner>, public ForcePrior
{
  public:
    explicit GeneralizedNewtonianViscousForce(BaseInnerRelation &inner_relation);
    virtual ~GeneralizedNewtonianViscousForce(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &mu_srd_;
};
using GeneralizedNewtonianViscousForceInner = GeneralizedNewtonianViscousForce<Inner<>>;

using BaseGeneralizedNewtonianViscousForceWithWall = InteractionWithWall<GeneralizedNewtonianViscousForce>;
template <>
class GeneralizedNewtonianViscousForce<Contact<Wall>> : BaseGeneralizedNewtonianViscousForceWithWall, public ForcePrior
{
  public:
    explicit GeneralizedNewtonianViscousForce(BaseContactRelation &contact_relation);
    virtual ~GeneralizedNewtonianViscousForce(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &mu_srd_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};
using GeneralizedNewtonianViscousForceWall = GeneralizedNewtonianViscousForce<Contact<Wall>>;
using GeneralizedNewtonianViscousForceWithWall = ComplexInteraction<GeneralizedNewtonianViscousForce<Inner<>, Contact<Wall>>>;

/**
 * @class VelocityGradient
 * @brief computes the shear rate of the fluid
 * @note this is needed for the ShearRateDependentViscosity (and therefore also the GeneralizedNewtonianViscousForce)
 */
template <typename... InteractionTypes>
class VelocityGradient;

template <class DataDelegationType>
class VelocityGradient<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit VelocityGradient(BaseRelationType &base_relation);
    virtual ~VelocityGradient(){};

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Matd> combined_velocity_gradient_;
};

template <>
class VelocityGradient<Inner<>> : public VelocityGradient<FluidDataInner>
{
  public:
    explicit VelocityGradient(BaseInnerRelation &inner_relation);
    virtual ~VelocityGradient(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &combined_velocity_gradient_;
};

template <>
class VelocityGradient<Contact<Wall>> : InteractionWithWall<VelocityGradient>
{
  public:
    explicit VelocityGradient(BaseContactRelation &contact_relation);
    virtual ~VelocityGradient(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
    StdLargeVec<Matd> &combined_velocity_gradient_;
};
using VelocityGradientWithWall = ComplexInteraction<VelocityGradient<Inner<>, Contact<Wall>>>;

/**
 * @class ShearRateDependentViscosity
 * @brief computes the viscosity based on the generalized Newtonian fluid model specified
 * @note needs VelocityGradient to work and is needed for GenearlizedNewtonianViscousForce
 */
class ShearRateDependentViscosity : public LocalDynamics, public FluidDataInner
{
  public:
    explicit ShearRateDependentViscosity(BaseInnerRelation &inner_relation);
    virtual ~ShearRateDependentViscosity(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &combined_velocity_gradient_;
    GeneralizedNewtonianFluid &generalized_newtonian_fluid_;
    StdLargeVec<Real> mu_srd_;
    StdLargeVec<Real> scalar_shear_rate_;
};

} // namespace fluid_dynamics
} // namespace SPH

#include "non_newtonian_dynamics.hpp"
#endif // NON_NEWTONIAN_DYNAMICS_H