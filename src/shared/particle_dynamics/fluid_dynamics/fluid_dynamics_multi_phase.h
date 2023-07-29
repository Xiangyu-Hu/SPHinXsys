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
 * @file fluid_dynamics_multi_phase.h
 * @brief Here, we define the algorithm classes for the dynamics involving multiple fluids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_MULTI_PHASE_H
#define FLUID_DYNAMICS_MULTI_PHASE_H

#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_complex.hpp"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateContact<BaseParticles, BaseParticles, DataDelegateEmptyBase>
    MultiPhaseContactData;
typedef DataDelegateContact<BaseParticles, BaseParticles> MultiPhaseData;
/**
 * @class ViscousAccelerationMultiPhase
 * @brief  the viscosity force induced acceleration
 */
class ViscousAccelerationMultiPhase : public ViscousAccelerationInner, public MultiPhaseContactData
{
  public:
    ViscousAccelerationMultiPhase(BaseInnerRelation &inner_relation,
                                  BaseContactRelation &contact_relation);
    explicit ViscousAccelerationMultiPhase(ComplexRelation &complex_relation);
    virtual ~ViscousAccelerationMultiPhase(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real> contact_mu_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};
using ViscousAccelerationMultiPhaseWithWall =
    BaseViscousAccelerationWithWall<ViscousAccelerationMultiPhase>;

/**
 * @class ViscousAccelerationMultiPhase
 * @brief Abstract base class for general multiphase fluid dynamics
 */
template <class RelaxationInnerType>
class RelaxationMultiPhase : public RelaxationInnerType, public MultiPhaseContactData
{
  public:
    RelaxationMultiPhase(BaseInnerRelation &inner_relation,
                         BaseContactRelation &contact_relation);
    virtual ~RelaxationMultiPhase(){};

  protected:
    StdVec<Fluid *> contact_fluids_;
    StdVec<StdLargeVec<Real> *> contact_p_, contact_rho_n_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};

/**
 * @class BaseMultiPhaseIntegration1stHalf
 * @brief  template class for multiphase pressure relaxation scheme
 */
template <class Integration1stHalfType>
class BaseMultiPhaseIntegration1stHalf : public RelaxationMultiPhase<Integration1stHalfType>
{
  public:
    BaseMultiPhaseIntegration1stHalf(BaseInnerRelation &inner_relation,
                                     BaseContactRelation &contact_relation);
    explicit BaseMultiPhaseIntegration1stHalf(ComplexRelation &complex_relation);
    virtual ~BaseMultiPhaseIntegration1stHalf(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    using CurrentRiemannSolver = decltype(Integration1stHalfType::riemann_solver_);
    StdVec<CurrentRiemannSolver> riemann_solvers_;
    virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
};
using MultiPhaseIntegration1stHalf = BaseMultiPhaseIntegration1stHalf<Integration1stHalf>;
using MultiPhaseIntegration1stHalfRiemann = BaseMultiPhaseIntegration1stHalf<Integration1stHalfRiemann>;

using MultiPhaseIntegration1stHalfWithWall =
    BaseIntegration1stHalfWithWall<MultiPhaseIntegration1stHalf>;
using MultiPhaseIntegration1stHalfRiemannWithWall =
    BaseIntegration1stHalfWithWall<MultiPhaseIntegration1stHalfRiemann>;
using ExtendMultiPhaseIntegration1stHalfRiemannWithWall =
    BaseExtendIntegration1stHalfWithWall<MultiPhaseIntegration1stHalfRiemann>;

/**
 * @class BaseMultiPhaseIntegration2ndHalf
 * @brief  template class pressure relaxation scheme with wall boundary
 */
template <class Integration2ndHalfType>
class BaseMultiPhaseIntegration2ndHalf : public RelaxationMultiPhase<Integration2ndHalfType>
{
  public:
    BaseMultiPhaseIntegration2ndHalf(BaseInnerRelation &inner_relation,
                                     BaseContactRelation &contact_relation);
    explicit BaseMultiPhaseIntegration2ndHalf(ComplexRelation &complex_relation);
    virtual ~BaseMultiPhaseIntegration2ndHalf(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    using CurrentRiemannSolver = decltype(Integration2ndHalfType::riemann_solver_);
    StdVec<CurrentRiemannSolver> riemann_solvers_;
};
using MultiPhaseIntegration2ndHalf = BaseMultiPhaseIntegration2ndHalf<Integration2ndHalf>;
using MultiPhaseIntegration2ndHalfRiemann = BaseMultiPhaseIntegration2ndHalf<Integration2ndHalfRiemann>;
using MultiPhaseIntegration2ndHalfWithWall = BaseMultiPhaseIntegration2ndHalf<MultiPhaseIntegration2ndHalf>;
using MultiPhaseIntegration2ndHalfRiemannWithWall = BaseIntegration2ndHalfWithWall<MultiPhaseIntegration2ndHalfRiemann>;

/**
 * @class MultiPhaseColorFunctionGradient
 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
 * TODO: Need a test cases for this.
 */
class MultiPhaseColorFunctionGradient : public LocalDynamics, public MultiPhaseData
{
  public:
    explicit MultiPhaseColorFunctionGradient(BaseContactRelation &contact_relation);
    virtual ~MultiPhaseColorFunctionGradient(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real rho0_;
    StdVec<Real> contact_rho0_;
    StdLargeVec<Real> &Vol_, &pos_div_;
    StdLargeVec<int> &surface_indicator_;
    StdLargeVec<Vecd> color_grad_, surface_norm_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_MULTI_PHASE_H