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
 * @file 	fluid_dynamics_complex.h
 * @brief 	Here, we define the algorithm classes for complex fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_COMPLEX_H
#define FLUID_DYNAMICS_COMPLEX_H

#include "base_fluid_dynamics.h"

#include "fluid_dynamics_inner.h"
#include "fluid_dynamics_inner.hpp"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateContact<BaseParticles, SolidParticles, DataDelegateEmptyBase> FluidWallData;
typedef DataDelegateContact<BaseParticles, SolidParticles> FSIContactData;
/**
 * @class InteractionWithWall
 * @brief Base class adding interaction with wall to general relaxation process
 */

template <template <class DataDelegationType> class BaseInteractionType>
class InteractionWithWall : public BaseInteractionType<FSIContactData>
{
  public:
    explicit InteractionWithWall(BaseContactRelation &wall_contact_relation);
    virtual ~InteractionWithWall(){};

  protected:
    StdVec<StdLargeVec<Vecd> *> wall_vel_ave_, wall_acc_ave_, wall_n_;
};

/**
 * @class BaseDensitySummationContact
 * @brief computing density by summation considering contribution from contact bodies
 */
class BaseDensitySummationContact : public BaseDensitySummation<FluidContactData>
{
  public:
    explicit BaseDensitySummationContact(BaseContactRelation &contact_relation);
    virtual ~BaseDensitySummationContact(){};

  protected:
    StdVec<Real> contact_inv_rho0_;
    StdVec<StdLargeVec<Real> *> contact_mass_;
    Real ContactSummation(size_t index_i);
};

/**
 * @class DensitySummationContact
 * @brief computing density by summation considering contribution from contact bodies
 */
class DensitySummationContact : public BaseDensitySummationContact
{
  public:
    explicit DensitySummationContact(BaseContactRelation &contact_relation)
        : BaseDensitySummationContact(contact_relation){};
    virtual ~DensitySummationContact(){};

    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class DensitySummationContactAdaptive
 * @brief computing density by summation considering  contribution from contact bodies
 */
class DensitySummationContactAdaptive : public BaseDensitySummationContact
{
  public:
    explicit DensitySummationContactAdaptive(BaseContactRelation &contact_relation);
    virtual ~DensitySummationContactAdaptive(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class ViscousWallBoundary
 * @brief TBD
 */
class ViscousWallBoundary : public InteractionWithWall<BaseViscousAcceleration>
{
  public:
    ViscousWallBoundary(BaseContactRelation &wall_contact_relation)
        : InteractionWithWall<BaseViscousAcceleration>(wall_contact_relation){};
    virtual ~ViscousWallBoundary(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class MomentumWallBoundary
 * @brief Wall boundary condition for solving the momentum equation
 */
template <class RiemannSolverType, class KernelCorrectionType>
class MomentumWallBoundary : public InteractionWithWall<BaseIntegration>
{
  public:
    MomentumWallBoundary(BaseContactRelation &wall_contact_relation);
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
    ExtendMomentumWallBoundary(BaseContactRelation &wall_contact_relation, Real penalty_strength = 1.0);
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
    ContinuityWallBoundary(BaseContactRelation &wall_contact_relation);
    virtual ~ContinuityWallBoundary(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using ContinuityWallBoundaryRiemann = ContinuityWallBoundary<AcousticRiemannSolver>;
using ContinuityWallBoundaryDissipativeRiemann = ContinuityWallBoundary<DissipativeRiemannSolver>;
/**
 * @class Oldroyd_BMomentumWallBoundary
 * @brief Dissipative Riemann solver is used here.
 */
class Oldroyd_BMomentumWallBoundary : public MomentumWallBoundaryDissipativeRiemann
{
  public:
    explicit Oldroyd_BMomentumWallBoundary(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BMomentumWallBoundary(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &tau_;
};

/**
 * @class Oldroyd_BContinuityWallBoundary
 * @brief Dissipative Riemann solver is used here too.
 */
class Oldroyd_BContinuityWallBoundary : public ContinuityWallBoundaryDissipativeRiemann
{
  public:
    explicit Oldroyd_BContinuityWallBoundary(BaseContactRelation &wall_contact_relation);
    virtual ~Oldroyd_BContinuityWallBoundary(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_COMPLEX_H
