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
 * @file multi_phase_dynamics.h
 * @brief Here, we define the algorithm classes for the dynamics involving multiple fluids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MULTI_PHASE_DYNAMICS_H
#define MULTI_PHASE_DYNAMICS_H

#include "fluid_integration.hpp"
#include "viscous_dynamics.hpp"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateContact<BaseParticles, BaseParticles> MultiPhaseContactData;

/**
 * @class MultiPhaseMomentumInterface
 * @brief Base class for general multiphase fluid dynamics
 */
template <class RiemannSolverType, class KernelCorrectionType>
class MultiPhaseMomentumInterface : public BaseIntegration<MultiPhaseContactData>
{
  public:
    MultiPhaseMomentumInterface(BaseContactRelation &contact_relation);
    virtual ~MultiPhaseMomentumInterface(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<KernelCorrectionType> correction_;
    StdVec<KernelCorrectionType> contact_corrections_;
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<StdLargeVec<Real> *> contact_p_;
};
using MultiPhaseMomentumInterfaceRiemann = MultiPhaseMomentumInterface<AcousticRiemannSolver, NoKernelCorrection>;

/**
 * @class MultiPhaseContinuityInterface
 * @brief  template class pressure relaxation scheme with wall boundary
 */
template <class RiemannSolverType>
class MultiPhaseContinuityInterface : public BaseIntegration<MultiPhaseContactData>
{
  public:
    MultiPhaseContinuityInterface(BaseContactRelation &contact_relation);
    virtual ~MultiPhaseContinuityInterface(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<StdLargeVec<Real> *> contact_p_;
    StdVec<StdLargeVec<Vecd> *> contact_vel_;
};
using MultiPhaseContinuityInterfaceRiemann = MultiPhaseContinuityInterface<AcousticRiemannSolver>;

/**
 * @class MultiPhaseColorFunctionGradient
 * @brief  indicate the particles near the interface of a fluid-fluid interaction and computing norm
 * TODO: Need a test cases for this.
 */
class MultiPhaseColorFunctionGradient : public LocalDynamics, public MultiPhaseContactData
{
  public:
    explicit MultiPhaseColorFunctionGradient(BaseContactRelation &contact_relation);
    virtual ~MultiPhaseColorFunctionGradient(){};

    void interaction(size_t index_i, Real dt = 0.0);

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
#endif // MULTI_PHASE_DYNAMICS_H