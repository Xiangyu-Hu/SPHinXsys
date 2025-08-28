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
 * @file 	continuum_integration.h
 * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
 * @details We consider here weakly compressible assumption to model elastic and
 * 			plastic materials with the updated Lagrangian framework.
 * @author	Shuaihao Zhang and Xiangyu Hu
 */
#ifndef CONTINUUM_INTEGRATION_H
#define CONTINUUM_INTEGRATION_H

#include "base_continuum_dynamics.h"
#include "constraint_dynamics.h"
#include "fluid_integration.hpp"
#include "general_continuum.h"
#include "general_continuum.hpp"
namespace SPH
{
namespace continuum_dynamics
{
class ContinuumInitialCondition : public LocalDynamics
{
  public:
    explicit ContinuumInitialCondition(SPHBody &sph_body);
    virtual ~ContinuumInitialCondition(){};

  protected:
    Vecd *pos_, *vel_;
    Mat3d *stress_tensor_3D_;
};

class AcousticTimeStep : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit AcousticTimeStep(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStep(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Fluid &fluid_;
    Real *rho_, *p_;
    Vecd *vel_;
    Real smoothing_length_min_;
    Real acousticCFL_;
};

template <class FluidDynamicsType>
class BaseIntegration1stHalf : public FluidDynamicsType
{
  public:
    explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalf(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *acc_shear_;
};
using Integration1stHalf = BaseIntegration1stHalf<fluid_dynamics::Integration1stHalfInnerNoRiemann>;
using Integration1stHalfRiemann = BaseIntegration1stHalf<fluid_dynamics::Integration1stHalfInnerRiemann>;

template <class DataDelegationType>
class BasePlasticIntegration : public fluid_dynamics::BaseIntegration<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit BasePlasticIntegration(BaseRelationType &base_relation);
    virtual ~BasePlasticIntegration(){};

  protected:
    PlasticContinuum &plastic_continuum_;
    Mat3d *stress_tensor_3D_, *strain_tensor_3D_, *stress_rate_3D_, *strain_rate_3D_;
    Matd *velocity_gradient_;
};

template <typename... InteractionTypes>
class PlasticIntegration1stHalf;

template <class RiemannSolverType>
class PlasticIntegration1stHalf<Inner<>, RiemannSolverType>
    : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit PlasticIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~PlasticIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
    virtual Vecd computeNonConservativeForce(size_t index_i);

  protected:
    RiemannSolverType riemann_solver_;
};
using PlasticIntegration1stHalfInnerNoRiemann = PlasticIntegration1stHalf<Inner<>, NoRiemannSolver>;
using PlasticIntegration1stHalfInnerRiemann = PlasticIntegration1stHalf<Inner<>, AcousticRiemannSolver>;

using BaseIntegrationWithWall = InteractionWithWall<BasePlasticIntegration>;

template <class RiemannSolverType>
class PlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit PlasticIntegration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~PlasticIntegration1stHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    virtual Vecd computeNonConservativeForce(size_t index_i);

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
using PlasticIntegration1stHalfWithWall = ComplexInteraction<PlasticIntegration1stHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
using PlasticIntegration1stHalfWithWallNoRiemann = PlasticIntegration1stHalfWithWall<NoRiemannSolver>;
using PlasticIntegration1stHalfWithWallRiemann = PlasticIntegration1stHalfWithWall<AcousticRiemannSolver>;

template <typename... InteractionTypes>
class PlasticIntegration2ndHalf;

template <class RiemannSolverType>
class PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>
    : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit PlasticIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~PlasticIntegration2ndHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    Real *Vol_, *mass_;
};
using PlasticIntegration2ndHalfInnerNoRiemann = PlasticIntegration2ndHalf<Inner<>, NoRiemannSolver>;
using PlasticIntegration2ndHalfInnerRiemann = PlasticIntegration2ndHalf<Inner<>, AcousticRiemannSolver>;

template <class RiemannSolverType>
class PlasticIntegration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit PlasticIntegration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~PlasticIntegration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
using PlasticIntegration2ndHalfWithWall = ComplexInteraction<PlasticIntegration2ndHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
using PlasticIntegration2ndHalfWithWallNoRiemann = PlasticIntegration2ndHalfWithWall<NoRiemannSolver>;
using PlasticIntegration2ndHalfWithWallRiemann = PlasticIntegration2ndHalfWithWall<AcousticRiemannSolver>;

class StressDiffusion : public BasePlasticIntegration<DataDelegateInner>
{
  public:
    explicit StressDiffusion(BaseInnerRelation &inner_relation);
    virtual ~StressDiffusion(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real zeta_ = 0.1, phi_; /*diffusion coefficient*/
    Real smoothing_length_, sound_speed_;
};

class ShearStressRelaxationHourglassControl1stHalf : public fluid_dynamics::BaseIntegration<DataDelegateInner>
{
  public:
    explicit ShearStressRelaxationHourglassControl1stHalf(BaseInnerRelation &inner_relation, Real xi = 4.0);
    virtual ~ShearStressRelaxationHourglassControl1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    GeneralContinuum &continuum_;
    Matd *shear_stress_, *velocity_gradient_, *strain_tensor_, *B_;
    Real *scale_penalty_force_, xi_;
};

class ShearStressRelaxationHourglassControl2ndHalf : public fluid_dynamics::BaseIntegration<DataDelegateInner>
{
  public:
    explicit ShearStressRelaxationHourglassControl2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~ShearStressRelaxationHourglassControl2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    GeneralContinuum &continuum_;
    Matd *shear_stress_, *velocity_gradient_;
    Vecd *acc_shear_, *acc_hourglass_;
    Real *scale_penalty_force_, G_;
};

class ShearStressRelaxationHourglassControl1stHalfJ2Plasticity : public ShearStressRelaxationHourglassControl1stHalf
{
  public:
    explicit ShearStressRelaxationHourglassControl1stHalfJ2Plasticity(BaseInnerRelation &inner_relation, Real xi = 0.2);
    virtual ~ShearStressRelaxationHourglassControl1stHalfJ2Plasticity(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    J2Plasticity &J2_plasticity_;
    Real *hardening_factor_;
};
} // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_INTEGRATION_H