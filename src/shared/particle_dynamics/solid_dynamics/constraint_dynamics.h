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

#include "all_particle_dynamics.h"
#include "all_simbody.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "general_constraint.h"
#include "general_reduce.h"
#include "solid_body.h"

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
    virtual ~SpringConstrain(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd stiffness_;
    StdLargeVec<Real> &mass_;
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
class PositionScaleSolidBody : public MotionConstraint<SPHBody>
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
class PositionTranslate : public MotionConstraint<DynamicsIdentifier>
{
  public:
    PositionTranslate(DynamicsIdentifier &identifier, Real start_time, Real end_time, Vecd translation)
        : MotionConstraint<DynamicsIdentifier>(identifier),
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
        return (this->pos0_[index_i] + translation_ - this->pos_[index_i]) * dt /
               (end_time_ - GlobalStaticVariables::physical_time_);
    };
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
    virtual void setupDynamics(Real dt = 0.0) override
    {
        velocity_correction_ =
            correction_matrix_ * compute_total_momentum_.exec(dt) / total_mass_;
    }

  public:
    explicit ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction = Vecd::Ones())
        : MotionConstraint<SPHBody>(sph_body), correction_matrix_(Matd::Identity()),
          compute_total_momentum_(sph_body, "Velocity")
    {
        for (int i = 0; i != Dimensions; ++i)
            correction_matrix_(i, i) = constrain_direction[i];
        ReduceDynamics<QuantitySummation<Real, SPHBody>> compute_total_mass_(sph_body, "Mass");
        total_mass_ = compute_total_mass_.exec();
    }
    virtual ~ConstrainSolidBodyMassCenter(){};

    void update(size_t index_i, Real dt = 0.0)
    {
        this->vel_[index_i] -= velocity_correction_;
    }
};

/**
 * @class ConstraintBySimBody
 * @brief Constrain by the motion computed from Simbody.
 */
template <class DynamicsIdentifier>
class ConstraintBySimBody : public MotionConstraint<DynamicsIdentifier>
{
  public:
    ConstraintBySimBody(DynamicsIdentifier &identifier,
                        SimTK::MultibodySystem &MBsystem,
                        SimTK::MobilizedBody &mobod,
                        SimTK::RungeKuttaMersonIntegrator &integ)
        : MotionConstraint<DynamicsIdentifier>(identifier),
          MBsystem_(MBsystem), mobod_(mobod), integ_(integ),
          n_(*this->particles_->template getVariableDataByName<Vecd>("NormalDirection")),
          n0_(*this->particles_->template registerSharedVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
          acc_(*this->particles_->template registerSharedVariable<Vecd>("Acceleration"))
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
         * base_particle_data_i.pos_  = (p_GB + r);
         */
        this->pos_[index_i] = degradeToVecd(SimTKToEigen(pos));
        this->vel_[index_i] = degradeToVecd(SimTKToEigen(vel));
        acc_[index_i] = degradeToVecd(SimTKToEigen(acc));

        SimTKVec3 n = (mobod_.getBodyRotation(*simbody_state_) * EigenToSimTK(upgradeToVec3d(n0_[index_i])));
        n_[index_i] = degradeToVecd(SimTKToEigen(n));
    };

  protected:
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    StdLargeVec<Vecd> &n_, &n0_, &acc_;
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
    : public BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>,
      public DataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &force_, &force_prior_, &pos_;
    SimTK::MultibodySystem &MBsystem_;
    SimTK::MobilizedBody &mobod_;
    SimTK::RungeKuttaMersonIntegrator &integ_;
    SimTKVec3 current_mobod_origin_location_;

  public:
    TotalForceForSimBody(DynamicsIdentifier &identifier,
                         SimTK::MultibodySystem &MBsystem,
                         SimTK::MobilizedBody &mobod,
                         SimTK::RungeKuttaMersonIntegrator &integ)
        : BaseLocalDynamicsReduce<ReduceSum<SimTK::SpatialVec>, DynamicsIdentifier>(identifier),
          DataDelegateSimple(identifier.getSPHBody()),
          force_(*particles_->registerSharedVariable<Vecd>("Force")),
          force_prior_(*particles_->getVariableDataByName<Vecd>("ForcePrior")),
          pos_(*particles_->getVariableDataByName<Vecd>("Position")),
          MBsystem_(MBsystem), mobod_(mobod), integ_(integ)
    {
        this->quantity_name_ = "TotalForceForSimBody";
    };

    virtual ~TotalForceForSimBody(){};

    virtual void setupDynamics(Real dt = 0.0) override
    {
        const SimTK::State *simbody_state = &integ_.getState();
        MBsystem_.realize(*simbody_state, SimTK::Stage::Acceleration);
        current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state);
    };

    SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd force = force_[index_i] + force_prior_[index_i];
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