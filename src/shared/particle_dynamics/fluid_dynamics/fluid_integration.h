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
 * @author	Xiangyu Hu
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
class FluidInitialCondition : public LocalDynamics
{
  public:
    explicit FluidInitialCondition(SPHBody &sph_body);
    virtual ~FluidInitialCondition(){};

  protected:
    Vecd *pos_, *vel_;
};

class ContinuumVolumeUpdate : public LocalDynamics
{
  public:
    explicit ContinuumVolumeUpdate(SPHBody &sph_body);
    virtual ~ContinuumVolumeUpdate(){};

    void update(size_t index_i, Real dt)
    {
        Vol_[index_i] = mass_[index_i] / rho_[index_i];
    }

  protected:
    Real *Vol_, *mass_, *rho_;
};

template <class DataDelegationType>
class BaseIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseIntegration(BaseRelationType &base_relation);
    virtual ~BaseIntegration(){};

  protected:
    Fluid &fluid_;
    Real *Vol_, *rho_, *mass_, *p_, *drho_dt_;
    Vecd *pos_, *vel_, *force_, *force_prior_;
};

template <typename... InteractionTypes>
class Integration1stHalf;

template <class RiemannSolverType, class KernelCorrectionType>
class Integration1stHalf<Inner<>, RiemannSolverType, KernelCorrectionType>
    : public BaseIntegration<DataDelegateInner>
{
  public:
    explicit Integration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Integration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};
using Integration1stHalfInnerNoRiemann = Integration1stHalf<Inner<>, NoRiemannSolver, NoKernelCorrection>;
using Integration1stHalfInnerRiemann = Integration1stHalf<Inner<>, AcousticRiemannSolver, NoKernelCorrection>;
using Integration1stHalfCorrectionInnerRiemann = Integration1stHalf<Inner<>, AcousticRiemannSolver, LinearGradientCorrection>;

// The following is used to avoid the C3200 error triggered in Visual Studio.
// Please refer: https://developercommunity.visualstudio.com/t/c-invalid-template-argument-for-template-parameter/831128
using BaseIntegrationWithWall = InteractionWithWall<BaseIntegration>;

template <class RiemannSolverType, class KernelCorrectionType>
class Integration1stHalf<Contact<Wall>, RiemannSolverType, KernelCorrectionType>
    : public BaseIntegrationWithWall
{
  public:
    explicit Integration1stHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Integration1stHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType, class KernelCorrectionType>
class Integration1stHalf<Contact<>, RiemannSolverType, KernelCorrectionType>
    : public BaseIntegration<DataDelegateContact>
{
  public:
    explicit Integration1stHalf(BaseContactRelation &contact_relation);
    virtual ~Integration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    StdVec<KernelCorrectionType> contact_corrections_;
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<Real *> contact_p_;
    StdVec<Real *> contact_Vol_;
};

template <class RiemannSolverType, class KernelCorrectionType>
using Integration1stHalfWithWall = ComplexInteraction<Integration1stHalf<Inner<>, Contact<Wall>>, RiemannSolverType, KernelCorrectionType>;

using Integration1stHalfWithWallNoRiemann = Integration1stHalfWithWall<NoRiemannSolver, NoKernelCorrection>;
using Integration1stHalfWithWallRiemann = Integration1stHalfWithWall<AcousticRiemannSolver, NoKernelCorrection>;
using Integration1stHalfCorrectionWithWallRiemann = Integration1stHalfWithWall<AcousticRiemannSolver, LinearGradientCorrection>;
using Integration1stHalfCorrectionForOpenBoundaryFlowWithWallRiemann = Integration1stHalfWithWall<AcousticRiemannSolver, LinearGradientCorrectionWithBulkScope>;

using MultiPhaseIntegration1stHalfWithWallRiemann =
    ComplexInteraction<Integration1stHalf<Inner<>, Contact<>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>;

template <typename... InteractionTypes>
class Integration2ndHalf;

template <class RiemannSolverType>
class Integration2ndHalf<Inner<>, RiemannSolverType>
    : public BaseIntegration<DataDelegateInner>
{
  public:
    typedef RiemannSolverType RiemannSolver;

    explicit Integration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Integration2ndHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    Real *mass_, *Vol_;
};
using Integration2ndHalfInnerRiemann = Integration2ndHalf<Inner<>, AcousticRiemannSolver>;
using Integration2ndHalfInnerNoRiemann = Integration2ndHalf<Inner<>, NoRiemannSolver>;
using Integration2ndHalfInnerDissipativeRiemann = Integration2ndHalf<Inner<>, DissipativeRiemannSolver>;

template <class RiemannSolverType>
class Integration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit Integration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~Integration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
class Integration2ndHalf<Contact<>, RiemannSolverType>
    : public BaseIntegration<DataDelegateContact>
{
  public:
    explicit Integration2ndHalf(BaseContactRelation &contact_relation);
    virtual ~Integration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdVec<RiemannSolverType> riemann_solvers_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_;
};

template <class RiemannSolverType>
using Integration2ndHalfWithWall = ComplexInteraction<Integration2ndHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;

using Integration2ndHalfWithWallNoRiemann = Integration2ndHalfWithWall<NoRiemannSolver>;
using Integration2ndHalfWithWallRiemann = Integration2ndHalfWithWall<AcousticRiemannSolver>;

using MultiPhaseIntegration2ndHalfWithWallRiemann =
    ComplexInteraction<Integration2ndHalf<Inner<>, Contact<>, Contact<Wall>>, AcousticRiemannSolver>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_INTEGRATION_H