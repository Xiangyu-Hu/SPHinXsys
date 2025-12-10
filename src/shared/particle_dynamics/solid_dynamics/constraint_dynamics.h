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
 * @file 	constraint_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONSTRAINT_DYNAMICS_H
#define CONSTRAINT_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "all_simbody.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "general_constraint.h"
#include "general_reduce.h"

namespace SPH
{
namespace solid_dynamics
{
/**@class SpringConstrain
 * @brief Constrain with a spring for each constrained particles to its original position.
 * //TODO: a test case is required for this class.
 */
class SpringConstrain : public MotionConstraint<BodyPartByParticle>
{
  public:
    SpringConstrain(BodyPartByParticle &body_part, Real stiffness);
    virtual ~SpringConstrain() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd stiffness_;
    Real *mass_;
    virtual Vecd getAcceleration(Vecd &disp, Real mass);
};

/**
 * @class 	PositionSolidBody
 * @brief 	Move a rigid body into a defined position in a given time interval,
 * 			can be considered as a quasi-static position driven boundary condition.
 * 			Note that, this constraint is not for a elastic solid body.
 */
class PositionSolidBody : public MotionConstraint<SPHBody>
{
  public:
    PositionSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Vecd pos_end_center);
    virtual ~PositionSolidBody() {};
    Vecd *GetParticlePos0() { return pos0_; };
    Vecd *GetParticlePosN() { return pos_; };
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real start_time_, end_time_;
    Real *physical_time_;
    Vecd pos_0_center_, pos_end_center_, translation_;
    Vecd getDisplacement(size_t index_i, Real dt);
};

/**
 * @class	PositionScaleSolidBody
 * @brief	Scale the body in a given time interval,
 * 			can be considered as a quasi-static position driven boundary condition.
 * 			Note that, this constraint is not for a elastic solid body.
 */
class PositionScaleSolidBody : public MotionConstraint<SPHBody>
{
  public:
    PositionScaleSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Real end_scale);
    virtual ~PositionScaleSolidBody() {};
    Vecd *GetParticlePos0() { return pos0_; };
    Vecd *GetParticlePosN() { return pos_; };
    virtual void update(size_t index_i, Real dt = 0.0);

  protected:
    Real start_time_, end_time_, end_scale_;
    Real *physical_time_;
    Vecd pos_0_center_;
    Vecd getDisplacement(size_t index_i, Real dt);
};

/**
 * @class	PositionTranslate
 * @brief	Translates the body in a given time interval
 * 			translation driven boundary condition; only moving the body; end position irrelevant;
 * 			Note that, this constraint is not for a elastic solid body.
 */
template <class DynamicsIdentifier>
class PositionTranslate : public MotionConstraint<DynamicsIdentifier>
{
  public:
    PositionTranslate(DynamicsIdentifier &identifier, Real start_time, Real end_time, Vecd translation);
    virtual ~PositionTranslate() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real start_time_, end_time_;
    Real *physical_time_;
    Vecd translation_;

    Vecd getDisplacement(size_t index_i, Real dt);
};
using TranslateSolidBody = PositionTranslate<SPHBody>;
using TranslateSolidBodyPart = PositionTranslate<BodyPartByParticle>;

/**
 * @class ConstrainSolidBodyMassCenter
 * @brief Constrain the mass center of a solid body.
 */
class ConstrainSolidBodyMassCenter : public MotionConstraint<SPHBody>
{
  private:
    Real total_mass_;
    Matd correction_matrix_;
    Vecd velocity_correction_;
    ReduceDynamics<QuantityMoment<Vecd, SPHBody>> compute_total_momentum_;

  protected:
    virtual void setupDynamics(Real dt = 0.0) override;

  public:
    explicit ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction = Vecd::Ones());
    virtual ~ConstrainSolidBodyMassCenter() {};
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ConstraintBySimBody
 * @brief Constrain by the motion computed from Simbody.
 */
template <class DynamicsIdentifier>
class ConstraintBySimBody : public MotionConstraint<DynamicsIdentifier>
{
  public:
    ConstraintBySimBody(DynamicsIdentifier &identifier, SimTK::MultibodySystem &MBsystem,
                        SimTK::MobilizedBody &mobod, SimTK::RungeKuttaMersonIntegrator &integ);
    virtual ~ConstraintBySimBody() {};
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
using ConstraintBodyBySimBody = ConstraintBySimBody<SPHBody>;
using ConstraintBodyPartBySimBody = ConstraintBySimBody<BodyPartByParticle>;

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

    virtual ~TotalForceForSimBody() {};

    virtual void setupDynamics(Real dt = 0.0) override;
    SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0);
};
using TotalForceOnBodyForSimBody = TotalForceForSimBody<SPHBody>;
using TotalForceOnBodyPartForSimBody = TotalForceForSimBody<BodyPartByParticle>;
} // namespace solid_dynamics
} // namespace SPH
#endif // CONSTRAINT_DYNAMICS_H