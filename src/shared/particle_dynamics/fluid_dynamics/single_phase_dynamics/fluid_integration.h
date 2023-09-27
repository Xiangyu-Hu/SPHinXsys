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
 * @file 	fluid_integration.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_INTEGRATION_H
#define FLUID_INTEGRATION_H

#include "base_fluid_dynamics.h"
#include "base_local_dynamics.h"
#include "riemann_solver.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class FluidInitialCondition
 * @brief  Set initial condition for a fluid body.
 * This is a abstract class to be override for case specific initial conditions
 */
class FluidInitialCondition : public LocalDynamics, public FluidDataSimple
{
  public:
    explicit FluidInitialCondition(SPHBody &sph_body);
    virtual ~FluidInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
};

/**
 * @class BaseIntegration
 * @brief Base class for all fluid relaxation schemes
 */
template <class DataDelegationType>
class BaseIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseIntegration(BaseRelationType &base_relation);
    virtual ~BaseIntegration(){};

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
};

/**
 * @class BaseIntegration1stHalfInner
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType, class KernelCorrectionType>
class BaseIntegration1stHalfInner : public BaseIntegration<FluidDataInner>
{
  public:
    explicit BaseIntegration1stHalfInner(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalfInner(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};
using Integration1stHalfInner = BaseIntegration1stHalfInner<NoRiemannSolver, NoKernelCorrection>;
using Integration1stHalfInnerRiemann = BaseIntegration1stHalfInner<AcousticRiemannSolver, NoKernelCorrection>;
using Integration1stHalfInnerDissipativeRiemann = BaseIntegration1stHalfInner<DissipativeRiemannSolver, NoKernelCorrection>;

/**
 * @class BaseIntegration2ndHalfInner
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalfInner : public BaseIntegration<FluidDataInner>
{
  public:
    explicit BaseIntegration2ndHalfInner(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration2ndHalfInner(){};
    void initialization(size_t index_i, Real dt = 0.0);
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Real> &Vol_, &mass_;
};
using Integration2ndHalfInner = BaseIntegration2ndHalfInner<NoRiemannSolver>;
/** define the mostly used density relaxation scheme using Riemann solver */
using Integration2ndHalfInnerRiemann = BaseIntegration2ndHalfInner<AcousticRiemannSolver>;
using Integration2ndHalfInnerDissipativeRiemann = BaseIntegration2ndHalfInner<DissipativeRiemannSolver>;

/**
 * @class MomentumWallBoundary
 * @brief Wall boundary condition for solving the momentum equation
 */
template <class RiemannSolverType, class KernelCorrectionType>
class MomentumWallBoundary : public InteractionWithWall<BaseIntegration>
{
  public:
    explicit MomentumWallBoundary(BaseContactRelation &wall_contact_relation);
    virtual ~MomentumWallBoundary(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

using MomentumWallBoundaryNoRiemann = MomentumWallBoundary<NoRiemannSolver, NoKernelCorrection>;
using MomentumWallBoundaryRiemann = MomentumWallBoundary<AcousticRiemannSolver, NoKernelCorrection>;
using MomentumWallBoundaryDissipativeRiemann = MomentumWallBoundary<DissipativeRiemannSolver, NoKernelCorrection>;
/**
 * @class ExtendMomentumWallBoundary
 * @brief Wall boundary conditions considering  wall penalty to prevent
 * particle penetration.
 */
template <class RiemannSolverType, class KernelCorrectionType>
class ExtendMomentumWallBoundary : public MomentumWallBoundary<RiemannSolverType, KernelCorrectionType>
{
  public:
    explicit ExtendMomentumWallBoundary(BaseContactRelation &wall_contact_relation, Real penalty_strength = 1.0);
    virtual ~ExtendMomentumWallBoundary(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real penalty_strength_;
};
using ExtendMomentumWallBoundaryRiemann = ExtendMomentumWallBoundary<AcousticRiemannSolver, NoKernelCorrection>;

/**
 * @class ContinuityWallBoundary
 * @brief Wall boundary condition for solving the continuity equation
 */
template <class RiemannSolverType>
class ContinuityWallBoundary : public InteractionWithWall<BaseIntegration>
{
  public:
    explicit ContinuityWallBoundary(BaseContactRelation &wall_contact_relation);
    virtual ~ContinuityWallBoundary(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using ContinuityWallBoundaryRiemann = ContinuityWallBoundary<AcousticRiemannSolver>;
using ContinuityWallBoundaryDissipativeRiemann = ContinuityWallBoundary<DissipativeRiemannSolver>;

class Integration1stHalfRiemannWithWall
    : public ComplexInteraction<Integration1stHalfInnerRiemann, MomentumWallBoundaryRiemann>
{
  public:
    explicit Integration1stHalfRiemannWithWall(ComplexRelation &fluid_wall_relation)
        : ComplexInteraction<Integration1stHalfInnerRiemann, MomentumWallBoundaryRiemann>(
              fluid_wall_relation.getInnerRelation(), fluid_wall_relation.getContactRelation()){};
};

class Integration2ndHalfRiemannWithWall
    : public ComplexInteraction<Integration2ndHalfInnerRiemann, ContinuityWallBoundaryRiemann>
{
  public:
    explicit Integration2ndHalfRiemannWithWall(ComplexRelation &fluid_wall_relation)
        : ComplexInteraction<Integration2ndHalfInnerRiemann, ContinuityWallBoundaryRiemann>(
              fluid_wall_relation.getInnerRelation(), fluid_wall_relation.getContactRelation()){};
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_INTEGRATION_H