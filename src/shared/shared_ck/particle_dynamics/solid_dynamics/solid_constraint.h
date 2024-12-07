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
    virtual ~ConstraintBySimBodyCK(){};

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
    SingularVariable<SimbodyState> sv_simbody_state_;

    void initializeSimbodyState(const SimTK::State &state);
    void updateSimbodyState(const SimTK::State &state);
};

template <class DynamicsIdentifier>
class ConstraintBySimBody : public MotionConstraint<DynamicsIdentifier>
{
  public:
    ConstraintBySimBody(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                        SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ);
    virtual ~ConstraintBySimBody(){};
    virtual void setupDynamics(Real dt = 0.0) override;
    void update(size_t index_i, Real dt = 0.0);

  protected:
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    SimbodyState simbody_state_;
    Vecd *n_, *n0_, *acc_;

    void initializeSimbodyState(const SimTK::State &state);
    void updateSimbodyState(const SimTK::State &state);
};

/**
 * @class TotalForceForSimBody
 * @brief Compute the force acting on the solid body part
 * for applying to simbody forces latter
 */
template <class DynamicsIdentifier>
class TotalForceForSimBody
    : public BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>
{
  protected:
    Vecd *force_, *force_prior_, *pos_;
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    SimTKVec3 current_mobod_origin_location_;

  public:
    TotalForceForSimBody(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                         SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ);

    virtual ~TotalForceForSimBody(){};

    virtual void setupDynamics(Real dt = 0.0) override;
    SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0);
};
using TotalForceOnBodyForSimBody = TotalForceForSimBody<SPHBody>;
using TotalForceOnBodyPartForSimBody = TotalForceForSimBody<BodyPartByParticle>;
} // namespace solid_dynamics
} // namespace SPH
#endif // SOLID_CONSTRAINT_H