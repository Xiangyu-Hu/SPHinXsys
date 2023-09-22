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
 * @file shear_dynamics.h
 * @brief shear dynamics in updated Lagrangian solid dynamics.
 * @author Shuaihao Zhang and Xiangyu Hu
 */

#ifndef SHEAR_DYNAMICS_H
#define SHEAR_DYNAMICS_H

#include "constraint_dynamics.h"
#include "continuum_particles.h"
#include "fluid_dynamics_complex.h"
#include "general_continuum.h"

namespace SPH
{
namespace continuum_dynamics
{
typedef DataDelegateSimple<ContinuumParticles> ContinuumDataSimple;
typedef DataDelegateInner<ContinuumParticles> ContinuumDataInner;

typedef DataDelegateSimple<PlasticContinuumParticles> PlasticContinuumDataSimple;
typedef DataDelegateInner<PlasticContinuumParticles> PlasticContinuumDataInner;

class ContinuumInitialCondition : public LocalDynamics, public PlasticContinuumDataSimple
{
  public:
    explicit ContinuumInitialCondition(SPHBody &sph_body);
    virtual ~ContinuumInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
    StdLargeVec<Mat3d> &stress_tensor_3D_;
};

class ContinuumAcousticTimeStepSize : public fluid_dynamics::AcousticTimeStepSize
{
  public:
    explicit ContinuumAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.5);
    virtual ~ContinuumAcousticTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
};

/**
 * @class BaseShearStressIntegration
 * @brief Evolution of shear stress
 */
template <class DataDelegationType>
class BaseShearStressIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseShearStressIntegration(BaseRelationType &base_relation);
    virtual ~BaseShearStressIntegration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<Matd> &velocity_gradient_, &shear_stress_;
    StdLargeVec<Real> &p_, &von_mises_stress_;
};

/**
 * @class ShearStressIntegration
 */
class ShearStressIntegration : public BaseShearStressIntegration<ContinuumDataInner>
{
  public:
    explicit ShearStressIntegration(BaseInnerRelation &inner_relation);
    virtual ~ShearStressIntegration(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    GeneralContinuum &continuum_;
};

/**
 * @class BaseShearStressAcceleration
 */
template <class DataDelegationType>
class BaseShearStressAcceleration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseShearStressAcceleration(BaseRelationType &base_relation);
    virtual ~BaseShearStressAcceleration(){};

  protected:
    StdLargeVec<Matd> &shear_stress_;
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &acc_prior_;
};

/**
 * @class ShearStressAcceleration
 */
class ShearStressAcceleration : public BaseShearStressAcceleration<ContinuumDataInner>
{
  public:
    explicit ShearStressAcceleration(BaseInnerRelation &inner_relation)
        : BaseShearStressAcceleration(inner_relation){};
    virtual ~ShearStressAcceleration(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShearStressAccelerationWithWall
 */
class ShearStressAccelerationWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressAcceleration>
{
  public:
    explicit ShearStressAccelerationWithWall(BaseContactRelation &wall_relation)
        : fluid_dynamics::InteractionWithWall<BaseShearStressAcceleration>(wall_relation){};
    virtual ~ShearStressAccelerationWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class PlasticShearStressIntegration
 */
class PlasticShearStressIntegration : public BaseShearStressIntegration<PlasticContinuumDataInner>
{
  public:
    explicit PlasticShearStressIntegration(BaseInnerRelation &inner_relation);
    virtual ~PlasticShearStressIntegration(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    StdLargeVec<Mat3d> &stress_tensor_3D_, &strain_tensor_3D_, &stress_rate_3D_, &strain_rate_3D_;
    StdLargeVec<Mat3d> &elastic_strain_tensor_3D_, &elastic_strain_rate_3D_;
    Real E_, nu_;
    StdLargeVec<Matd> &shear_stress_;
};

/**
 * @class ShearAccelerationIntegration
 * @brief This designed for hourglass free formulation
 */
class ShearAccelerationIntegration : public LocalDynamics, public ContinuumDataInner
{
  public:
    explicit ShearAccelerationIntegration(BaseInnerRelation &inner_relation);
    virtual ~ShearAccelerationIntegration(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    GeneralContinuum &continuum_;
    Real G_, smoothing_length_;
    StdLargeVec<Vecd> &vel_, &acc_prior_, acc_shear_;
    StdLargeVec<Real> &rho_;
};

/**
 * @class StressDiffusion
 */
class StressDiffusion : public LocalDynamics, public PlasticContinuumDataInner
{
  public:
    explicit StressDiffusion(BaseInnerRelation &inner_relation, SharedPtr<Gravity> gravity_ptr, int axis);
    virtual ~StressDiffusion(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    PlasticContinuum &plastic_continuum_;
    int axis_;
    Real zeta_ = 0.1, rho0_, gravity_, smoothing_length_;
    Real phi_, diffusion_coeff_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Mat3d> &stress_tensor_3D_, &stress_rate_3D_;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // SHEAR_DYNAMICS_H
