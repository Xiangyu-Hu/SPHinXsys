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
 * @file 	constraint_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef CONSTRAINT_DYNAMICS_H
#define CONSTRAINT_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "general_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "all_simbody.h"

namespace SPH
{
namespace solid_dynamics
{
//----------------------------------------------------------------------
//		for general solid dynamics
//----------------------------------------------------------------------
typedef DataDelegateSimple<SolidParticles> SolidDataSimple;
typedef DataDelegateInner<SolidParticles> SolidDataInner;

/**
 * @class	BaseMotionConstraint
 * @brief	Base class for constraining with prescribed motion.
 * 			Exact motion function will be defined in derive class.
 * 			Note that, we do not impose acceleration, so that this constraint
 * 			must be imposed after updating particle velocity by forces
 * 			and before updating particle position.
 * 			TODO: to clarify the treatment of particle position,
 * 			how to achieve consistency between velocity and position constraints.
 */
template <class DynamicsIdentifier>
class BaseMotionConstraint : public BaseLocalDynamics<DynamicsIdentifier>, public SolidDataSimple
{
  public:
    explicit BaseMotionConstraint(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier), SolidDataSimple(identifier.getSPHBody()),
          pos_(particles_->pos_), pos0_(particles_->pos0_),
          n_(particles_->n_), n0_(particles_->n0_),
          vel_(particles_->vel_), acc_(particles_->acc_){};

    virtual ~BaseMotionConstraint(){};

  protected:
    StdLargeVec<Vecd> &pos_, &pos0_;
    StdLargeVec<Vecd> &n_, &n0_;
    StdLargeVec<Vecd> &vel_, &acc_;
};

/**@class FixConstraint
 * @brief Constraint with zero velocity.
 */
template <class DynamicsIdentifier>
class FixConstraint : public BaseMotionConstraint<DynamicsIdentifier>
{
  public:
    explicit FixConstraint(DynamicsIdentifier &identifier)
        : BaseMotionConstraint<DynamicsIdentifier>(identifier){};
    virtual ~FixConstraint(){};

    void update(size_t index_i, Real dt = 0.0) { this->vel_[index_i] = Vecd::Zero(); };
};
using FixBodyConstraint = FixConstraint<SPHBody>;
using FixBodyPartConstraint = FixConstraint<BodyPartByParticle>;

/**@class SpringConstrain
 * @brief Constrain with a spring for each constrained particles to its original position.
 * //TODO: a test case is required for this class.
 */
class SpringConstrain : public BaseMotionConstraint<BodyPartByParticle>
{
  public:
    SpringConstrain(BodyPartByParticle &body_part, Real stiffness);
    virtual ~SpringConstrain(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &mass_;
    Vecd stiffness_;
    virtual Vecd getAcceleration(Vecd &disp, Real mass);
};

/**
 * @class 	PositionSolidBody
 * @brief 	Move a rigid body into a defined position in a given time interval,
 * 			can be considered as a quasi-static position driven boundary condition.
 * 			Note that, this constraint is not for a elastic solid body.
 */
class PositionSolidBody : public BaseMotionConstraint<SPHBody>
{
  public:
    PositionSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Vecd pos_end_center);
    virtual ~PositionSolidBody(){};
    StdLargeVec<Vecd> &GetParticlePos0() { return pos0_; };
    StdLargeVec<Vecd> &GetParticlePosN() { return pos_; };
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real start_time_, end_time_;
    Vecd pos_0_center_, pos_end_center_, translation_;
    Vecd getDisplacement(size_t index_i, Real dt);
};

/**
 * @class	PositionScaleSolidBody
 * @brief	Scale the body in a given time interval,
 * 			can be considered as a quasi-static position driven boundary condition.
 * 			Note that, this constraint is not for a elastic solid body.
 */
class PositionScaleSolidBody : public BaseMotionConstraint<SPHBody>
{
  public:
    PositionScaleSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Real end_scale);
    virtual ~PositionScaleSolidBody(){};
    StdLargeVec<Vecd> &GetParticlePos0() { return pos0_; };
    StdLargeVec<Vecd> &GetParticlePosN() { return pos_; };
    virtual void update(size_t index_i, Real dt = 0.0);

  protected:
    Real start_time_, end_time_, end_scale_;
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
class PositionTranslate : public BaseMotionConstraint<DynamicsIdentifier>
{
  public:
    PositionTranslate(DynamicsIdentifier &identifier, Real start_time, Real end_time, Vecd translation)
        : BaseMotionConstraint<DynamicsIdentifier>(identifier),
          start_time_(start_time), end_time_(end_time), translation_(translation){};
    virtual ~PositionTranslate(){};
    void update(size_t index_i, Real dt = 0.0)
    {
        // only apply in the defined time period
        if (GlobalStaticVariables::physical_time_ >= start_time_ && GlobalStaticVariables::physical_time_ <= end_time_)
        {
            // displacement from the initial position, 0.5x because it's executed twice
            this->pos_[index_i] += 0.5 * getDisplacement(index_i, dt);
            this->vel_[index_i] = Vecd::Zero();
        }
    };

  protected:
    Real start_time_, end_time_;
    Vecd translation_;
    Vecd getDisplacement(size_t index_i, Real dt)
    {
        Vecd displacement = Vecd::Zero();
        displacement = (this->pos0_[index_i] + translation_ - this->pos_[index_i]) * dt /
                       (end_time_ - GlobalStaticVariables::physical_time_);
        return displacement;
    };
};
using TranslateSolidBody = PositionTranslate<SPHBody>;
using TranslateSolidBodyPart = PositionTranslate<BodyPartByParticle>;

/**
 * @class FixedInAxisDirection
 * @brief Constrain the velocity of a solid body part.
 */
class FixedInAxisDirection : public BaseMotionConstraint<BodyPartByParticle>
{
  public:
    FixedInAxisDirection(BodyPartByParticle &body_part, Vecd constrained_axises = Vecd::Zero());
    virtual ~FixedInAxisDirection(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd constrain_matrix_;
};

/**
 * @class ConstrainSolidBodyMassCenter
 * @brief Constrain the mass center of a solid body.
 */
class ConstrainSolidBodyMassCenter : public LocalDynamics, public SolidDataSimple
{
  private:
    Real total_mass_;
    Matd correction_matrix_;
    Vecd velocity_correction_;
    StdLargeVec<Vecd> &vel_;
    ReduceDynamics<QuantityMoment<Vecd>> compute_total_momentum_;

  protected:
    virtual void setupDynamics(Real dt = 0.0) override;

  public:
    explicit ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction = Vecd::Ones());
    virtual ~ConstrainSolidBodyMassCenter(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ConstraintBySimBody
 * @brief Constrain by the motion computed from Simbody.
 */
template <class DynamicsIdentifier>
class ConstraintBySimBody : public BaseMotionConstraint<DynamicsIdentifier>
{
  public:
    ConstraintBySimBody(DynamicsIdentifier &identifier,
                        SimTK::MultibodySystem &MBsystem,
                        SimTK::MobilizedBody &mobod,
                        SimTK::RungeKuttaMersonIntegrator &integ)
        : BaseMotionConstraint<DynamicsIdentifier>(identifier),
          MBsystem_(MBsystem), mobod_(mobod), integ_(integ)
    {
        simbody_state_ = &integ_.getState();
        MBsystem_.realize(*simbody_state_, SimTK::Stage::Acceleration);
        initial_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
    };
    virtual ~ConstraintBySimBody(){};

    virtual void setupDynamics(Real dt = 0.0) override
    {
        simbody_state_ = &integ_.getState();
        MBsystem_.realize(*simbody_state_, SimTK::Stage::Acceleration);
    };
    void update(size_t index_i, Real dt = 0.0)
    {
        /** Change to SimTK::Vector. */
        SimTKVec3 rr, pos, vel, acc;
        rr = EigenToSimTK(upgradeToVec3d(this->pos0_[index_i])) - initial_mobod_origin_location_;
        mobod_.findStationLocationVelocityAndAccelerationInGround(*simbody_state_, rr, pos, vel, acc);
        /** this is how we calculate the particle position in after transform of MBbody.
         * const SimTK::Rotation&  R_GB = mobod_.getBodyRotation(simbody_state);
         * const SimTKVec3&      p_GB = mobod_.getBodyOriginLocation(simbody_state);
         * const SimTKVec3 r = R_GB * rr; // re-express station vector p_BS in G (15 flops)
         * base_particle_data_i.pos_ = (p_GB + r);
         */
        degradeToVecd(SimTKToEigen(pos), this->pos_[index_i]);
        degradeToVecd(SimTKToEigen(vel), this->vel_[index_i]);

        SimTKVec3 n = (mobod_.getBodyRotation(*simbody_state_) * EigenToSimTK(upgradeToVec3d(this->n0_[index_i])));
        degradeToVecd(SimTKToEigen(n), this->n_[index_i]);
    };

  protected:
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    const SimTK::State *simbody_state_;
    SimTKVec3 initial_mobod_origin_location_;
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
    : public BaseLocalDynamicsReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>,
      public SolidDataSimple
{
  protected:
    StdLargeVec<Real> &mass_;
    StdLargeVec<Vecd> &acc_, &acc_prior_, &pos_;
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    const SimTK::State *simbody_state_;
    SimTKVec3 current_mobod_origin_location_;

  public:
    TotalForceForSimBody(DynamicsIdentifier &identifier,
                         SimTK::MultibodySystem &MBsystem,
                         SimTK::MobilizedBody &mobod,
                         SimTK::RungeKuttaMersonIntegrator &integ)
        : BaseLocalDynamicsReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>(
              identifier, SimTK::SpatialVec(SimTKVec3(0), SimTKVec3(0))),
          SolidDataSimple(identifier.getSPHBody()), mass_(particles_->mass_),
          acc_(particles_->acc_), acc_prior_(particles_->acc_prior_),
          pos_(particles_->pos_),
          MBsystem_(MBsystem), mobod_(mobod), integ_(integ)
    {
        this->quantity_name_ = "TotalForceForSimBody";
    };

    virtual ~TotalForceForSimBody(){};

    virtual void setupDynamics(Real dt = 0.0) override
    {
        simbody_state_ = &integ_.getState();
        MBsystem_.realize(*simbody_state_, SimTK::Stage::Acceleration);
        current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
    };

    SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd force = (acc_[index_i] + acc_prior_[index_i]) * mass_[index_i];
        SimTKVec3 force_from_particle = EigenToSimTK(upgradeToVec3d(force));
        SimTKVec3 displacement = EigenToSimTK(upgradeToVec3d(pos_[index_i])) - current_mobod_origin_location_;
        SimTKVec3 torque_from_particle = SimTK::cross(displacement, force_from_particle);

        return SimTK::SpatialVec(torque_from_particle, force_from_particle);
    };
};
using TotalForceOnBodyForSimBody = TotalForceForSimBody<SPHBody>;
using TotalForceOnBodyPartForSimBody = TotalForceForSimBody<BodyPartByParticle>;
} // namespace solid_dynamics
} // namespace SPH
#endif // CONSTRAINT_DYNAMICS_H