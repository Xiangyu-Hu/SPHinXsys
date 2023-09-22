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
typedef DataDelegateContact<BaseParticles, BaseParticles> MultiPhaseContactData;
/**
 * @class ViscousAccelerationMultiPhase
 * @brief  the viscosity force induced acceleration
 */
class ViscousAccelerationMultiPhase : public BaseViscousAcceleration<MultiPhaseContactData>
{
  public:
    explicit ViscousAccelerationMultiPhase(BaseContactRelation &contact_relation);
    virtual ~ViscousAccelerationMultiPhase(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<Real> contact_mu_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};

/**
 * @class MultiPhaseMomentumInterface
 * @brief Abstract base class for general multiphase fluid dynamics
 */
template <class RiemannSolverType, class InterfacePressureType>
class MultiPhaseMomentumInterface : public BaseIntegration<MultiPhaseContactData>
{
  public:
    MultiPhaseMomentumInterface(BaseContactRelation &contact_relation);
    virtual ~MultiPhaseMomentumInterface(){};

  protected:
    InterfacePressureType interface_pressure_;
    StdVec<Fluid *> contact_fluids_;
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<StdLargeVec<Real> *> contact_p_, contact_rho_n_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};
using MultiPhaseIntegration1stHalfRiemann = MultiPhaseMomentumInterface<AcousticRiemannSolver, PlainInterfacePressure>;

using MultiPhaseIntegration1stHalfWithWall =
    MomentumWallBoundary<MultiPhaseIntegration1stHalf>;
using MultiPhaseIntegration1stHalfRiemannWithWall =
    MomentumWallBoundary<MultiPhaseIntegration1stHalfRiemann>;
using ExtendMultiPhaseIntegration1stHalfRiemannWithWall =
    ExtendMomentumWallBoundary<MultiPhaseIntegration1stHalfRiemann>;

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
using MultiPhaseIntegration2ndHalfRiemannWithWall = ContinuityWallBoundary<MultiPhaseIntegration2ndHalfRiemann>;

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
    StdLargeVec<int> &indicator_;
    StdLargeVec<Vecd> color_grad_, surface_norm_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_MULTI_PHASE_H