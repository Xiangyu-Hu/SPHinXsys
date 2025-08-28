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
 * @file 	eulerian_fluid_integration.h
 * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_FLUID_INTEGRATION_H
#define EULERIAN_FLUID_INTEGRATION_H

#include "fluid_integration.hpp"
#include "riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
template <class DataDelegationType>
class EulerianIntegration : public BaseIntegration<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit EulerianIntegration(BaseRelationType &base_relation);
    virtual ~EulerianIntegration() {};

  protected:
    Vecd *mom_, *dmom_dt_;
    Real *dmass_dt_, *Vol_;
};

template <typename... InteractionTypes>
class EulerianIntegration1stHalf;

template <class RiemannSolverType>
class EulerianIntegration1stHalf<Inner<>, RiemannSolverType>
    : public EulerianIntegration<DataDelegateInner>
{
  public:
    explicit EulerianIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration1stHalf(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration1stHalf(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration1stHalf() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration1stHalfInnerRiemann = EulerianIntegration1stHalf<Inner<>, AcousticRiemannSolver>;

using BaseEulerianIntegrationWithWall = InteractionWithWall<EulerianIntegration>;
template <class RiemannSolverType>
class EulerianIntegration1stHalf<Contact<Wall>, RiemannSolverType>
    : public BaseEulerianIntegrationWithWall
{
  public:
    EulerianIntegration1stHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration1stHalf(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration1stHalf(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration1stHalf() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration1stHalfWithWallRiemann =
    ComplexInteraction<EulerianIntegration1stHalf<Inner<>, Contact<Wall>>, AcousticRiemannSolver>;

template <typename... InteractionTypes>
class EulerianIntegration2ndHalf;

template <class RiemannSolverType>
class EulerianIntegration2ndHalf<Inner<>, RiemannSolverType>
    : public EulerianIntegration<DataDelegateInner>
{
  public:
    typedef RiemannSolverType RiemannSolver;

    explicit EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration2ndHalf(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration2ndHalf(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration2ndHalf() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfInnerRiemann = EulerianIntegration2ndHalf<Inner<>, AcousticRiemannSolver>;

/**
 * @class EulerianIntegration2ndHalfWithWall
 * @brief template density relaxation scheme with using  Riemann solver.
 */
template <class RiemannSolverType>
class EulerianIntegration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseEulerianIntegrationWithWall
{
  public:
    EulerianIntegration2ndHalf(BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0);
    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianIntegration2ndHalf(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianIntegration2ndHalf(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~EulerianIntegration2ndHalf() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfWithWallRiemann =
    ComplexInteraction<EulerianIntegration2ndHalf<Inner<>, Contact<Wall>>, AcousticRiemannSolver>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_INTEGRATION_H
