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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	solid_constraint.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SOLID_CONSTRAINT_H
#define SOLID_CONSTRAINT_H

#include "all_simbody.h"
#include "general_constraint_ck.h"
#include "general_reduce_ck.h"
#include "solid_body.h"

namespace SPH
{
namespace solid_dynamics
{

template <class DynamicsIdentifier>
class ConstraintBySimBodyCK : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    explicit ConstraintBySimBodyCK(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                                   SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ);
    virtual ~ConstraintBySimBodyCK() {};
    virtual void setupDynamics(Real dt = 0.0) override;
    SingularVariable<SimbodyState> *svSimbodyState() { return sv_simbody_state_; };

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *pos_, *pos0_, *vel_;
        Vecd *n_, *n0_, *acc_;
        SimbodyState *simbody_state_;
    };

  protected:
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    DiscreteVariable<Vecd> *dv_pos_, *dv_pos0_, *dv_vel_;
    DiscreteVariable<Vecd> *dv_n_, *dv_n0_, *dv_acc_;
    SingularVariable<SimbodyState> *sv_simbody_state_;
    SimTKVec3 sim_tk_initial_origin_location_;
};
using ConstraintBodyBySimBodyCK = ConstraintBySimBodyCK<SPHBody>;
using ConstraintBodyPartBySimBodyCK = ConstraintBySimBodyCK<BodyPartByParticle>;

template <class DynamicsIdentifier>
class TotalForceForSimBodyCK
    : public BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>
{

  public:
    TotalForceForSimBodyCK(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                           SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ);

    virtual ~TotalForceForSimBodyCK() {};
    virtual void setupDynamics(Real dt = 0.0) override;

    class ReduceKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *force_, *force_prior_, *pos_;
        Vec3d *current_origin_location_;
    };

  protected:
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    DiscreteVariable<Vecd> *dv_force_, *dv_force_prior_, *dv_pos_;
    SingularVariable<Vec3d> *sv_current_origin_location_;
};
using TotalForceOnBodyForSimBodyCK = TotalForceForSimBodyCK<SPHBody>;
using TotalForceOnBodyPartForSimBodyCK = TotalForceForSimBodyCK<BodyPartByParticle>;
} // namespace solid_dynamics
} // namespace SPH
#endif // SOLID_CONSTRAINT_H