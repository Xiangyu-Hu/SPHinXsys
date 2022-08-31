/**
 * @file 	constraint_dynamics.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "constraint_dynamics.h"

#include <numeric>

using namespace SimTK;

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		BaseMotionConstraint::BaseMotionConstraint(SPHBody &sph_body)
			: LocalDynamics(sph_body), SolidDataSimple(sph_body_),
			  pos_(particles_->pos_), pos0_(particles_->pos0_),
			  n_(particles_->n_), n0_(particles_->n0_),
			  vel_(particles_->vel_), acc_(particles_->acc_) {}
		//=================================================================================================//
		BaseMotionConstraint::BaseMotionConstraint(BodyPartByParticle &body_part)
			: BaseMotionConstraint(body_part.getSPHBody()) {}
		//=================================================================================================//
		SpringConstrain::SpringConstrain(BodyPartByParticle &body_part, Real stiffness)
			: BaseMotionConstraint(body_part), mass_(particles_->mass_), stiffness_(stiffness) {}
		//=================================================================================================//
		Vecd SpringConstrain::getAcceleration(Vecd &disp, Real mass)
		{
			Vecd spring_force(0);
			for (int i = 0; i < disp.size(); i++)
			{
				spring_force[i] = -stiffness_[i] * disp[i] / mass;
			}
			return spring_force;
		}
		//=================================================================================================//
		void SpringConstrain::update(size_t index_i, Real dt)
		{
			Vecd displacement = pos_[index_i] - pos0_[index_i];
			vel_[index_i] += dt * getAcceleration(displacement, mass_[index_i]);
		}
		//=================================================================================================//
		PositionSolidBody::
			PositionSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Vecd pos_end_center)
			: BaseMotionConstraint(sph_body),
			  start_time_(start_time), end_time_(end_time), pos_end_center_(pos_end_center)
		{
			BoundingBox bounds = sph_body.getBodyShapeBounds();
			pos_0_center_ = (bounds.first + bounds.second) * 0.5;
			translation_ = pos_end_center_ - pos_0_center_;
		}
		//=================================================================================================//
		Vecd PositionSolidBody::getDisplacement(size_t index_i, Real dt)
		{
			Vecd displacement;
			// displacement from the initial position
			Vecd pos_final = pos0_[index_i] + translation_;
			displacement = (pos_final - pos_[index_i]) * dt /
						   (end_time_ - GlobalStaticVariables::physical_time_);

			return displacement;
		}
		//=================================================================================================//
		void PositionSolidBody::update(size_t index_i, Real dt)
		{
			// only apply in the defined time period
			if (GlobalStaticVariables::physical_time_ >= start_time_ &&
				GlobalStaticVariables::physical_time_ <= end_time_)
			{
				pos_[index_i] = pos_[index_i] + getDisplacement(index_i, dt); // displacement from the initial position
				vel_[index_i] = Vecd(0);
			}
		}
		//=================================================================================================//
		PositionScaleSolidBody::
			PositionScaleSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Real end_scale)
			: BaseMotionConstraint(sph_body),
			  start_time_(start_time), end_time_(end_time), end_scale_(end_scale)
		{
			BoundingBox bounds = sph_body.getBodyShapeBounds();
			pos_0_center_ = (bounds.first + bounds.second) * 0.5;
		}
		//=================================================================================================//
		Vecd PositionScaleSolidBody::getDisplacement(size_t index_i, Real dt)
		{
			Vecd displacement(0);
			// displacement from the initial position
			Vecd pos_final = pos_0_center_ + end_scale_ * (pos0_[index_i] - pos_0_center_);
			displacement = (pos_final - pos_[index_i]) * dt /
						   (end_time_ - GlobalStaticVariables::physical_time_);
			return displacement;
		}
		//=================================================================================================//
		void PositionScaleSolidBody::update(size_t index_i, Real dt)
		{
			// only apply in the defined time period
			if (GlobalStaticVariables::physical_time_ >= start_time_ &&
				GlobalStaticVariables::physical_time_ <= end_time_)
			{
				pos_[index_i] = pos_[index_i] + getDisplacement(index_i, dt); // displacement from the initial position
				vel_[index_i] = Vecd(0);
			}
		}
		//=================================================================================================//
		TranslateSolidBody::
			TranslateSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Vecd translation)
			: BaseMotionConstraint(sph_body),
			  start_time_(start_time), end_time_(end_time), translation_(translation) {}
		//=================================================================================================//
		TranslateSolidBody::
			TranslateSolidBody(BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd translation)
			: TranslateSolidBody(body_part.getSPHBody(), start_time, end_time, translation){};
		//=================================================================================================//
		Vecd TranslateSolidBody::getDisplacement(size_t index_i, Real dt)
		{
			Vecd displacement(0);
			displacement = (pos0_[index_i] + translation_ - pos_[index_i]) * dt / (end_time_ - GlobalStaticVariables::physical_time_);
			return displacement;
		}
		//=================================================================================================//
		void TranslateSolidBody::update(size_t index_i, Real dt)
		{
			// only apply in the defined time period
			if (GlobalStaticVariables::physical_time_ >= start_time_ && GlobalStaticVariables::physical_time_ <= end_time_)
			{
				pos_[index_i] = pos_[index_i] + 0.5 * getDisplacement(index_i, dt); // displacement from the initial position, 0.5x because it's executed twice
				vel_[index_i] = Vecd(0);
			}
		}
		//=================================================================================================//
		ConstrainSolidBodyMassCenter::
			ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction)
			: LocalDynamics(sph_body), SolidDataSimple(sph_body),
			  correction_matrix_(Matd(1.0)), vel_(particles_->vel_),
			  compute_total_momentum_(sph_body, "Velocity")
		{
			for (int i = 0; i != Dimensions; ++i)
				correction_matrix_[i][i] = constrain_direction[i];
			ReduceDynamics<QuantitySummation<Real>> compute_total_mass_(sph_body, "MassiveMeasure");
			total_mass_ = compute_total_mass_.parallel_exec();
		}
		//=================================================================================================//
		void ConstrainSolidBodyMassCenter::setupDynamics(Real dt)
		{
			velocity_correction_ =
				correction_matrix_ * compute_total_momentum_.parallel_exec(dt) / total_mass_;
		}
		//=================================================================================================//
		void ConstrainSolidBodyMassCenter::update(size_t index_i, Real dt)
		{
			vel_[index_i] -= velocity_correction_;
		}
		//=================================================================================================//
		ConstraintBySimBody::ConstraintBySimBody(SPHBody &sph_body,
												 SimTK::MultibodySystem &MBsystem,
												 SimTK::MobilizedBody &mobod,
												 SimTK::Force::DiscreteForces &force_on_bodies,
												 SimTK::RungeKuttaMersonIntegrator &integ)
			: BaseMotionConstraint(sph_body),
			  MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			initial_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
		ConstraintBySimBody::ConstraintBySimBody(BodyPartByParticle &body_part,
												 SimTK::MultibodySystem &MBsystem,
												 SimTK::MobilizedBody &mobod,
												 SimTK::Force::DiscreteForces &force_on_bodies,
												 SimTK::RungeKuttaMersonIntegrator &integ)
			: ConstraintBySimBody(body_part.getSPHBody(), MBsystem, mobod, force_on_bodies, integ){};
		//=================================================================================================//
		void ConstraintBySimBody::setupDynamics(Real dt)
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
		}
		//=================================================================================================//
		TotalForceForSimBody::TotalForceForSimBody(SPHBody &sph_body,
												   SimTK::MultibodySystem &MBsystem,
												   SimTK::MobilizedBody &mobod,
												   SimTK::Force::DiscreteForces &force_on_bodies,
												   SimTK::RungeKuttaMersonIntegrator &integ)
			: LocalDynamicsReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>(sph_body, SpatialVec(Vec3(0), Vec3(0))),
			  SolidDataSimple(sph_body), mass_(particles_->mass_),
			  acc_(particles_->acc_), acc_prior_(particles_->acc_prior_),
			  pos_(particles_->pos_),
			  MBsystem_(MBsystem), mobod_(mobod), force_on_bodies_(force_on_bodies), integ_(integ)
		{
			quantity_name_ = "TotalForceForSimBody";
		}
		//=================================================================================================//
		TotalForceForSimBody::TotalForceForSimBody(BodyPartByParticle &body_part,
												   SimTK::MultibodySystem &MBsystem,
												   SimTK::MobilizedBody &mobod,
												   SimTK::Force::DiscreteForces &force_on_bodies,
												   SimTK::RungeKuttaMersonIntegrator &integ)
			: TotalForceForSimBody(body_part.getSPHBody(), MBsystem, mobod, force_on_bodies, integ) {}
		//=================================================================================================//
		void TotalForceForSimBody::setupDynamics(Real dt)
		{
			simbody_state_ = &integ_.getState();
			MBsystem_.realize(*simbody_state_, Stage::Acceleration);
			current_mobod_origin_location_ = mobod_.getBodyOriginLocation(*simbody_state_);
		}
		//=================================================================================================//
	}
}
