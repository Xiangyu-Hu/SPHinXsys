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
// ********************** //
// Herschel Bulkley Fluid
// ********************** //
//!
// class ShearRateDependentViscosity : public PairGeomAverageFixed<Real>
// {
//   public:
//     ShearRateDependentViscosity(Real visc1, Real visc2)
//         : PairGeomAverageFixed<Real>(
//               visc1,
//               visc2){};
//     explicit ShearRateDependentViscosity(Real visc)
//         : ShearRateDependentViscosity(visc, visc){};
//     virtual ~ShearRateDependentViscosity(){};

// };

template <typename... InteractionTypes>
class HerschelBulkleyAcceleration;

template <class DataDelegationType>
class HerschelBulkleyAcceleration<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit HerschelBulkleyAcceleration(BaseRelationType &base_relation);
    virtual ~HerschelBulkleyAcceleration(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_;
    StdLargeVec<Vecd> &vel_, &viscous_force_;
    HerschelBulkleyFluid &herschel_bulkley_fluid_;
    Real smoothing_length_;
};

template <>
class HerschelBulkleyAcceleration<Inner<>> : public HerschelBulkleyAcceleration<FluidDataInner>, public ForcePrior
{
  public:
    explicit HerschelBulkleyAcceleration(BaseInnerRelation &inner_relation);
    virtual ~HerschelBulkleyAcceleration(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &mu_srd_;
};
using HerschelBulkleyAccelerationInner = HerschelBulkleyAcceleration<Inner<>>;

using BaseHerschelBulkleyAccelerationWithWall = InteractionWithWall<HerschelBulkleyAcceleration>;
template <>
class HerschelBulkleyAcceleration<Contact<Wall>> : BaseHerschelBulkleyAccelerationWithWall, public ForcePrior
{
  public:
    explicit HerschelBulkleyAcceleration(BaseContactRelation &contact_relation);
    virtual ~HerschelBulkleyAcceleration(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &mu_srd_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};
using HerschelBulkleyAccelerationWall = HerschelBulkleyAcceleration<Contact<Wall>>;
using ViscousShearRateDependent = ComplexInteraction<HerschelBulkleyAcceleration<Inner<>, Contact<Wall>>>;

/**
 * @class VelocityGradientInner
 * @brief  computes the shear rate of the fluid
 */

class VelocityGradientInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit VelocityGradientInner(BaseInnerRelation &inner_relation);
    virtual ~VelocityGradientInner(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Matd> combined_velocity_gradient_;
};

/**
 * @class VelocityGradientContact
 * @brief  computes the shear rate of the fluid at the boundary
 */

class VelocityGradientContact : public LocalDynamics, public FluidContactData
{
  public:
    explicit VelocityGradientContact(BaseContactRelation &contact_relation);
    virtual ~VelocityGradientContact(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Matd> &combined_velocity_gradient_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};

class ShearRateDependentViscosity : public LocalDynamics, public FluidDataInner
{
  public:
    explicit ShearRateDependentViscosity(BaseInnerRelation &inner_relation);
    virtual ~ShearRateDependentViscosity(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &combined_velocity_gradient_;
    HerschelBulkleyFluid &herschel_bulkley_fluid_;
    StdLargeVec<Real> mu_srd_;
    StdLargeVec<Real> scalar_shear_rate_;
};

} // namespace fluid_dynamics
} // namespace SPH

#include "non_newtonian_dynamics.hpp"
#endif // NON_NEWTONIAN_DYNAMICS_H