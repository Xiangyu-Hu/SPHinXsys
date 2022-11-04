/* -----------------------------------------------------------------------------*
 *                               SPHinXsys                                      *
 * -----------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle    *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for       *
 * physical accurate simulation and aims to model coupled industrial dynamic    *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH      *
 * (smoothed particle hydrodynamics), a meshless computational method using     *
 * particle discretization.                                                     *
 *                                                                              *
 * SPHinXsys is partially funded by German Research Foundation                  *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,               *
 * HU1527/12-1 and HU1527/12-4.                                                 *
 *                                                                              *
 * Portions copyright (c) 2017-2022 Technical University of Munich and          *
 * the authors' affiliations.                                                   *
 *                                                                              *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may      *
 * not use this file except in compliance with the License. You may obtain a    *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.           *
 *                                                                              *
 * -----------------------------------------------------------------------------*/
/**
 * @file 	constraint_dynamics.h
 * @brief 	Here, we define the algorithm classes for solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#ifndef CONSTRAINT_DYNAMICS_H
#define CONSTRAINT_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"

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
		class BaseMotionConstraint : public LocalDynamics, public SolidDataSimple
		{
		public:
			explicit BaseMotionConstraint(SPHBody &sph_body);
			explicit BaseMotionConstraint(BodyPartByParticle &body_part);

			virtual ~BaseMotionConstraint(){};

		protected:
			StdLargeVec<Vecd> &pos_, &pos0_;
			StdLargeVec<Vecd> &n_, &n0_;
			StdLargeVec<Vecd> &vel_, &acc_;
		};

		/**@class FixConstraint
		 * @brief Constraint with zero velocity.
		 */
		class FixConstraint : public BaseMotionConstraint
		{
		public:
			explicit FixConstraint(SPHBody &sph_body) : BaseMotionConstraint(sph_body){};
			explicit FixConstraint(BodyPartByParticle &body_part) : BaseMotionConstraint(body_part){};
			virtual ~FixConstraint(){};

			void update(size_t index_i, Real dt = 0.0) { vel_[index_i] = Vecd(0); };
		};

		/**@class SpringConstrain
		 * @brief Constrain with a spring for each constrained particles to its original position.
		 * //TODO: a test case is required for this class.
		 */
		class SpringConstrain : public BaseMotionConstraint
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
		class PositionSolidBody : public BaseMotionConstraint
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
		class PositionScaleSolidBody : public BaseMotionConstraint
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
		 * @class	TranslateSolidBody
		 * @brief	Translates the body in a given time interval
		 * 			translation driven boundary condition; only moving the body; end position irrelevant;
		 * 			Note that, this constraint is not for a elastic solid body.
		 */
		class TranslateSolidBody : public BaseMotionConstraint
		{
		public:
			TranslateSolidBody(SPHBody &sph_body, Real start_time, Real end_time, Vecd translation);
			TranslateSolidBody(BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd translation);
			virtual ~TranslateSolidBody(){};
			void update(size_t index_i, Real dt = 0.0);

		protected:
			Real start_time_, end_time_;
			Vecd translation_;
			Vecd getDisplacement(size_t index_i, Real dt);
		};

		/**
		 * @class FixedInAxisDirection
		 * @brief Constrain the velocity of a solid body part.
		 */
		class FixedInAxisDirection : public BaseMotionConstraint
		{
		public:
			FixedInAxisDirection(BodyPartByParticle &body_part, Vecd constrained_axises = Vecd(0))
				: BaseMotionConstraint(body_part), constrain_matrix_(Matd(1.0))
			{
				for (int k = 0; k != Dimensions; ++k)
					constrain_matrix_[k][k] = constrained_axises[k];
			};
			virtual ~FixedInAxisDirection(){};

			void update(size_t index_i, Real dt = 0.0)
			{
				vel_[index_i] = constrain_matrix_ * vel_[index_i];
			};

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
			explicit ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction = Vecd(1.0));
			virtual ~ConstrainSolidBodyMassCenter(){};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class ConstraintBySimBody
		 * @brief Constrain by the motion computed from Simbody.
		 */
		class ConstraintBySimBody : public BaseMotionConstraint
		{
		public:
			ConstraintBySimBody(SPHBody &sph_body,
								SimTK::MultibodySystem &MBsystem,
								SimTK::MobilizedBody &mobod,
								SimTK::Force::DiscreteForces &force_on_bodies,
								SimTK::RungeKuttaMersonIntegrator &integ);
			ConstraintBySimBody(BodyPartByParticle &body_part,
								SimTK::MultibodySystem &MBsystem,
								SimTK::MobilizedBody &mobod,
								SimTK::Force::DiscreteForces &force_on_bodies,
								SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~ConstraintBySimBody(){};

			virtual void setupDynamics(Real dt = 0.0) override;
			void virtual update(size_t index_i, Real dt = 0.0);

		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody &mobod_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			const SimTK::State *simbody_state_;
			Vec3d initial_mobod_origin_location_;
		};

		/**
		 * @class TotalForceForSimBody
		 * @brief Compute the force acting on the solid body part
		 * for applying to simbody forces latter
		 */
		class TotalForceForSimBody
			: public LocalDynamicsReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>,
			  public SolidDataSimple
		{
		protected:
			StdLargeVec<Real> &mass_;
			StdLargeVec<Vecd> &acc_, &acc_prior_, &pos_;
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody &mobod_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			const SimTK::State *simbody_state_;
			Vec3d current_mobod_origin_location_;

		public:
			TotalForceForSimBody(SPHBody &sph_body,
								 SimTK::MultibodySystem &MBsystem,
								 SimTK::MobilizedBody &mobod,
								 SimTK::Force::DiscreteForces &force_on_bodies,
								 SimTK::RungeKuttaMersonIntegrator &integ);
			TotalForceForSimBody(BodyPartByParticle &body_part,
								 SimTK::MultibodySystem &MBsystem,
								 SimTK::MobilizedBody &mobod,
								 SimTK::Force::DiscreteForces &force_on_bodies,
								 SimTK::RungeKuttaMersonIntegrator &integ);

			virtual ~TotalForceForSimBody(){};

			virtual void setupDynamics(Real dt = 0.0) override;
			SimTK::SpatialVec reduce(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif // CONSTRAINT_DYNAMICS_H