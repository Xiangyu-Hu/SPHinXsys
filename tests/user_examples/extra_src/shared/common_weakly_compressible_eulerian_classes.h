/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2023 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	common_weakly_compressible_eulerian_classes.h
 * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef COMMON_WEAKLY_COMPRESSIBLE_EULERIAN_CLASSES_H
#define COMMON_WEAKLY_COMPRESSIBLE_EULERIAN_CLASSES_H

#include "compressible_fluid.h"
#include "fluid_body.h"
#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
//----------------------------------------------------------------------
//	Riemann Solver classes.
//----------------------------------------------------------------------
/**
 * @struct NoRiemannSolverInWCEulerianMethod
 * @brief  NO RiemannSolver for weakly-compressible flow in Eulerian method for weakly-compressible flow.
 */
class NoRiemannSolverInWCEulerianMethod
{
    Fluid &fluid_i_, &fluid_j_;

  public:
    NoRiemannSolverInWCEulerianMethod(Fluid &fluid_i, Fluid &fluid_j)
        : fluid_i_(fluid_i), fluid_j_(fluid_j){};
    FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij);
};

/**
 * @struct AcousticRiemannSolverInEulerianMethod
 * @brief  Acoustic RiemannSolver for weakly-compressible flow in Eulerian method.
 */
class AcousticRiemannSolverInEulerianMethod
{
    Fluid &fluid_i_, &fluid_j_;
    Real limiter_parameter_;

  public:
    AcousticRiemannSolverInEulerianMethod(Fluid &compressible_fluid_i, Fluid &compressible_fluid_j, Real limiter_parameter = 15.0)
        : fluid_i_(compressible_fluid_i), fluid_j_(compressible_fluid_j), limiter_parameter_(limiter_parameter){};
    FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij);
};

/**
 * @class BaseIntegration1stHalf
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType>
class BaseIntegration1stHalf : public fluid_dynamics::BaseIntegration
{
  public:
    explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    virtual ~BaseIntegration1stHalf(){};
    Real limiter_input_;
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &acc_prior_;
    StdLargeVec<Vecd> mom_, dmom_dt_;
};
/** define the mostly used pressure relaxation scheme using Riemann solver */
using Integration1stHalfAcousticRiemann = BaseIntegration1stHalf<AcousticRiemannSolverInEulerianMethod>;

/**
 * @class BaseIntegration1stHalfWithWall
 * @brief  template class pressure relaxation scheme with wall boundary
 */
template <class BaseIntegration1stHalfType>
class BaseIntegration1stHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseIntegration1stHalfType>
{
  public:
    // template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    BaseIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation, Real limiter_parameter = 15.0)
        : fluid_dynamics::InteractionWithWall<BaseIntegration1stHalfType>(wall_contact_relation, base_body_relation), limiter_input_(limiter_parameter){};
    explicit BaseIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation)
        : BaseIntegration1stHalfWithWall(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~BaseIntegration1stHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
    Real &limiter_input_;
};
using Integration1stHalfAcousticRiemannWithWall = BaseIntegration1stHalfWithWall<Integration1stHalfAcousticRiemann>;

/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalf : public fluid_dynamics::BaseIntegration
{
  public:
    explicit BaseIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    virtual ~BaseIntegration2ndHalf(){};
    Real limiter_input_;
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using Integration2ndHalfAcousticRiemann = BaseIntegration2ndHalf<AcousticRiemannSolverInEulerianMethod>;

/**
 * @class BaseIntegration2ndHalfWithWall
 * @brief template density relaxation scheme with using  Riemann solver.
 */
template <class BaseIntegration2ndHalfType>
class BaseIntegration2ndHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseIntegration2ndHalfType>
{
  public:
    // template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    BaseIntegration2ndHalfWithWall(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation, Real limiter_parameter = 15.0)
        : fluid_dynamics::InteractionWithWall<BaseIntegration2ndHalfType>(wall_contact_relation, base_body_relation), limiter_input_(limiter_parameter){};
    explicit BaseIntegration2ndHalfWithWall(ComplexRelation &fluid_wall_relation)
        : BaseIntegration2ndHalfWithWall(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~BaseIntegration2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
    Real &limiter_input_;
};
using Integration2ndHalfAcousticRiemannWithWall = BaseIntegration2ndHalfWithWall<Integration2ndHalfAcousticRiemann>;

//----------------------------------------------------------------------
//	Non-Reflective Boundary
//----------------------------------------------------------------------
class NonReflectiveBoundaryVariableCorrection : public LocalDynamics, public DataDelegateInner<BaseParticles>
{
  public:
    NonReflectiveBoundaryVariableCorrection(BaseInnerRelation &inner_relation);
    virtual ~NonReflectiveBoundaryVariableCorrection(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
    Real rho_farfield_, sound_speed_;
    Vecd vel_farfield_;
    StdLargeVec<Real> &rho_, &p_, &Vol_;
    StdLargeVec<Vecd> &vel_, &mom_, &pos_;
    StdLargeVec<Vecd> n_;
    StdLargeVec<Real> inner_weight_summation_, rho_average_, vel_normal_average_;
    StdLargeVec<Vecd> vel_tangential_average_, vel_average_;
    StdLargeVec<int> &indicator_;
    StdLargeVec<int> surface_inner_particle_indicator_;
};
} // namespace SPH
#endif // COMMON_WEAKLY_COMPRESSIBLE_EULERIAN_CLASSES_H