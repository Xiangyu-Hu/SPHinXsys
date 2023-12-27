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
 * @file 	gerenal_eulerian_fluid_dynamics.h
 * @brief 	Here, we define the general eulerian fluid dynamcis classes for fluid dynamics.
 * @author	Bo Zhang and Xiangyu Hu
 */
#ifndef GENERAL_EULERIAN_FLUID_DYNAMICS_H
#define GENERAL_EULERIAN_FLUID_DYNAMICS_H

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
  * @struct DissipationState
  * @brief struct for stored dissipation state of Riemann solver for eulerian fluid dynamics.
  */
struct DissipationState
{
    Matd momentum_dissipation_;
    Vecd density_dissipation_;

    DissipationState(Matd momentum_dissipation = Matd::Zero(), Vecd density_dissipation = Vecd::Zero()) : 
        momentum_dissipation_(momentum_dissipation), density_dissipation_(density_dissipation){};
};

/**
  * @struct NoRiemannSolverGE
  * @brief
  */
class NoRiemannSolverGE
{
  public:
    NoRiemannSolverGE(Fluid &fluid_i, Fluid &fluid_j) : fluid_i_(fluid_i), fluid_j_(fluid_j){};
    Vec2d getBoundingWaveSpeeds(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij);
    DissipationState getDissipationState(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij);

  protected:
    Fluid &fluid_i_, &fluid_j_;
};

/**
  * @ struct HLLERiemannSolver
  * @ brief
  */
class HLLERiemannSolver : public NoRiemannSolverGE
{
  public:
    HLLERiemannSolver(Fluid &fluid_i, Fluid &fluid_j) : NoRiemannSolverGE(fluid_i, fluid_j){};
    Vec2d getBoundingWaveSpeeds(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij);
    DissipationState getDissipationState(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij);
};

/**
  * @class ICEIntegration1stHalf
  * @brief Template class for pressure relaxation scheme with the Riemann solver
  * as template variable.
  */
template <class RiemannSolverType>
class ICEIntegration1stHalf : public BaseIntegration
{
  public:
    explicit ICEIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~ICEIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Vecd> &acc_prior_;
    StdLargeVec<Vecd> mom_, dmom_dt_;
};
/** define the mostly used pressure relaxation scheme using Riemann solver */
using ICEIntegration1stHalfNoRiemann = ICEIntegration1stHalf<NoRiemannSolverGE>;
using ICEIntegration1stHalfHLLERiemann = ICEIntegration1stHalf<HLLERiemannSolver>;

/**
  * @class ICEIntegration1stHalfWithWall
  * @brief Template class pressure relaxation schemem with wall boundary.
  */
template <class EulerianIntegration1stHalfType>
class ICEIntegration1stHalfWithWall : public InteractionWithWall<EulerianIntegration1stHalfType>
{
    public:
    // template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    ICEIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation)
        : InteractionWithWall<EulerianIntegration1stHalfType>(wall_contact_relation, base_body_relation){};
    explicit ICEIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation)
        : ICEIntegration1stHalfWithWall(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~ICEIntegration1stHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using ICEIntegration1stHalfNoRiemannWithWall = ICEIntegration1stHalfWithWall<ICEIntegration1stHalfNoRiemann>;
using ICEIntegration1stHalfHLLERiemannWithWall = ICEIntegration1stHalfWithWall<ICEIntegration1stHalfHLLERiemann>;

/**
  * @class EulerianIntegration2ndHalf
  * @brief Template denstiy relaxation scheme with different Riemann solver.
  */
template <class RiemannSolverType>
class ICEIntegration2ndHalf : public BaseIntegration
{
  public:
    explicit ICEIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~ICEIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

    protected:
    RiemannSolverType riemann_solver_;
};
using ICEIntegration2ndHalfNoRiemann = ICEIntegration2ndHalf<NoRiemannSolverGE>;
using ICEIntegration2ndHalfHLLERiemann = ICEIntegration2ndHalf<HLLERiemannSolver>;

/**
  * @class ICEIntegration2ndHalfWithWall
  * @brief Template density relaxation scheme with using Riemann solver.
  */
template <class EulerianIntegration2ndHalfType>
class ICEIntegration2ndHalfWithWall : public InteractionWithWall<EulerianIntegration2ndHalfType>
{
  public:
    // template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    ICEIntegration2ndHalfWithWall(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation)
        : InteractionWithWall<EulerianIntegration2ndHalfType>(wall_contact_relation, base_body_relation){};
    explicit ICEIntegration2ndHalfWithWall(ComplexRelation &fluid_wall_relation)
        : ICEIntegration2ndHalfWithWall(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~ICEIntegration2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using ICEIntegration2ndHalfNoRiemannWithWall = ICEIntegration2ndHalfWithWall<ICEIntegration2ndHalfNoRiemann>;
using ICEIntegration2ndHalfHLLERiemannWithWall = ICEIntegration2ndHalfWithWall<ICEIntegration2ndHalfHLLERiemann>;

/**
  * @class CheckTaylorGreenVortexFlow
  * @brief Check the numerical result of taylor green vortex and gives the error.
  */
class CheckTaylorGreenVortexFlow : public LocalDynamics, public FluidDataInner
{
public:
    explicit CheckTaylorGreenVortexFlow(BaseInnerRelation& inner_relation);
    virtual ~CheckTaylorGreenVortexFlow() {};

    void update(size_t index_i, Real dt = 0.0);

protected:
    StdLargeVec<Vecd>& pos_, & vel_;
    StdLargeVec<Real> velocity_error_norm_;
};

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
#endif // GENERAL_EULERIAN_FLUID_DYNAMICS_H