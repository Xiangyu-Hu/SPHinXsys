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
 * @file 	eulerian_fluid_dynamics.h
 * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_FLUID_DYNAMICS_H
#define EULERIAN_FLUID_DYNAMICS_H

#include "compressible_fluid.h"
#include "fluid_body.h"
#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @struct EulerianAcousticRiemannSolver
 * @brief  Acoustic RiemannSolver for Eulerian weakly-compressible flow.
 */
class EulerianAcousticRiemannSolver
{
    Fluid &fluid_i_, &fluid_j_;
    Real limiter_parameter_;

  public:
    EulerianAcousticRiemannSolver(Fluid &fluid_i, Fluid &fluid_j, Real limiter_parameter = 15.0)
        : fluid_i_(fluid_i), fluid_j_(fluid_j), limiter_parameter_(limiter_parameter){};
    FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij);
};

/**
 * @class WeaklyCompressibleFluidInitialCondition
 * @brief  Set initial condition for a Eulerian weakly compressible fluid body.
 * This is a abstract class to be override for case specific initial conditions
 */
class WeaklyCompressibleFluidInitialCondition : public FluidInitialCondition
{
  public:
    explicit WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), rho_(particles_->rho_),
          p_(*particles_->getVariableByName<Real>("Pressure")),
          mom_(*particles_->getVariableByName<Vecd>("Momentum")){};

  protected:
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &mom_;
};

/**
 * @class EulerianIntegration1stHalfInner
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType>
class EulerianIntegration1stHalfInner : public BaseIntegration<FluidDataInner>
{
  public:
    explicit EulerianIntegration1stHalfInner(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    virtual ~EulerianIntegration1stHalfInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Vecd> mom_, dmom_dt_;
};
/** define the mostly used pressure relaxation scheme using Riemann solver */
using EulerianIntegration1stHalfInnerAcousticRiemann = EulerianIntegration1stHalfInner<EulerianAcousticRiemannSolver>;

/**
 * @class EulerianIntegration1stHalfWithWall
 * @brief  template class pressure relaxation scheme with wall boundary
 */
template <class RiemannSolverType>
class EulerianIntegration1stHalfWithWall : public InteractionWithWall<BaseIntegration>
{
  public:
    EulerianIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0);
    virtual ~EulerianIntegration1stHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Vecd> &dmom_dt_;
};
using EulerianIntegration1stHalfWithWallAcousticRiemann = EulerianIntegration1stHalfWithWall<EulerianAcousticRiemannSolver>;

/**
 * @class EulerianIntegration2ndHalf
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class EulerianIntegration2ndHalf : public BaseIntegration<FluidDataInner>
{
  public:
    explicit EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    virtual ~EulerianIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfAcousticRiemann = EulerianIntegration2ndHalf<EulerianAcousticRiemannSolver>;

/**
 * @class EulerianIntegration2ndHalfWithWall
 * @brief template density relaxation scheme with using  Riemann solver.
 */
template <class RiemannSolverType>
class EulerianIntegration2ndHalfWithWall : public InteractionWithWall<BaseIntegration>
{
  public:
    EulerianIntegration2ndHalfWithWall(BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0);
    virtual ~EulerianIntegration2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfWithWallAcousticRiemann = EulerianIntegration2ndHalfWithWall<EulerianAcousticRiemannSolver>;

/**
 * @class SmearedSurfaceIndication
 * @brief Indication of the particles which are within cut-off radius of surface particles.
 */
class SmearedSurfaceIndication : public LocalDynamics, public FluidDataInner
{
  public:
    explicit SmearedSurfaceIndication(BaseInnerRelation &inner_relation);
    virtual ~SmearedSurfaceIndication(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<int> &indicator_;
    StdLargeVec<int> &smeared_surface_;
};

/**
 * @class NonReflectiveBoundaryCorrection
 * @brief Implement Eulerian non-reflective boundary condition at free surface particles.
 */
class NonReflectiveBoundaryCorrection : public LocalDynamics, public DataDelegateInner<BaseParticles>
{
  public:
    NonReflectiveBoundaryCorrection(BaseInnerRelation &inner_relation);
    virtual ~NonReflectiveBoundaryCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
    Real rho_farfield_, sound_speed_;
    Vecd vel_farfield_;
    StdLargeVec<Real> &rho_, &p_, &Vol_;
    StdLargeVec<Vecd> &vel_, &mom_, &pos_;
    StdLargeVec<Real> inner_weight_summation_, rho_average_, vel_normal_average_;
    StdLargeVec<Vecd> vel_tangential_average_, vel_average_;
    StdLargeVec<int> &indicator_, smeared_surface_;
    StdLargeVec<Vecd> &n_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_DYNAMICS_H