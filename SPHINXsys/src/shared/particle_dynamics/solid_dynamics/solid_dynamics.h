/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	solid_dynamics.h
* @brief 	Here, we define the algorithm classes for solid dynamics. 
* @details 	We consider here a weakly compressible solids.   
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once
#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "weakly_compressible_fluid.h"
#include "base_kernel.h"
#include "fluid_dynamics.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for general solid dynamics 
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidBody, SolidParticles, Solid> SolidDataDelegateSimple;
		typedef DataDelegateInner<SolidBody, SolidParticles, Solid> SolidDataDelegateInner;
		typedef DataDelegateContact<SolidBody, SolidParticles, Solid, FluidBody, FluidParticles, Fluid> FSIDataDelegateContact;
		typedef DataDelegateContact<SolidBody, SolidParticles, Solid, SolidBody, SolidParticles, Solid> ContactDynamicsDataDelegate;

		/**
		 * @class SolidDynamicsInitialCondition
		 * @brief  set initial condition for solid fluid body
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class SolidDynamicsInitialCondition :
			public ParticleDynamicsSimple, public SolidDataDelegateSimple
		{
		public:
			SolidDynamicsInitialCondition(SolidBody* body) :
				ParticleDynamicsSimple(body), SolidDataDelegateSimple(body) {};
			virtual ~SolidDynamicsInitialCondition() {};
		};

		/**
		* @class SummationContactDensity
		* @brief Computing the summation density due to solid-solid contact model.
		*/
		class SummationContactDensity :
			public PartInteractionDynamicsByParticle, public ContactDynamicsDataDelegate
		{
		public:
			SummationContactDensity(SolidBodyContactRelation* soild_body_contact_relation);
			virtual ~SummationContactDensity() {};
		protected:
			StdLargeVec<Real>& mass_, & contact_density_;
			StdVec<StdLargeVec<Real>*> contact_mass_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ContactForceAcceleration
		* @brief Computing the contact force.
		*/
		class ContactForce :
			public PartInteractionDynamicsByParticle, public ContactDynamicsDataDelegate
		{
		public:
			ContactForce(SolidBodyContactRelation* soild_body_contact_relation);
			virtual ~ContactForce() {};
		protected:
			StdLargeVec<Real>& contact_density_, & Vol_, & mass_;
			StdLargeVec<Vecd>& dvel_dt_others_, & contact_force_;
			StdVec<StdLargeVec<Real>*> contact_contact_density_, contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ContactForceFromFriction
		* @brief  Computing the contact-force contribution from friction.
		*/
		class ContactForceFromFriction :
			public InteractionDynamicsSplitting, public ContactDynamicsDataDelegate
		{
		public:
			ContactForceFromFriction(SPHBodyContactRelation* body_contact_relation,
				StdLargeVec<Vecd>& vel_n, Real eta);
			virtual ~ContactForceFromFriction() {};
		protected:
			StdLargeVec<Real>& Vol_, & mass_;
			StdLargeVec<Vecd>& contact_force_, & vel_n_;
			StdVec<StdLargeVec<Real>*> contact_Vol_, contact_mass_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_n_, contact_contact_force_;
			Real eta_; /**< friction coefficient */

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FluidViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidViscousForceOnSolid :
			public InteractionDynamics, public FSIDataDelegateContact
		{
		public:
			FluidViscousForceOnSolid(SPHBodyContactRelation* body_contact_relation);
			virtual ~FluidViscousForceOnSolid() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& vel_ave_, & viscous_force_from_fluid_;
			StdVec<StdLargeVec<Real>*> contact_Vol_, contact_rho_n_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_n_;
			StdVec<Real> mu_;
			StdVec<Real> smoothing_length_;

			//dynamics of a particle
			//to be realized in specific algorithms
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FluidAngularConservativeViscousForceOnSolid
		* @brief Computing the viscous force from the fluid
		*/
		class FluidAngularConservativeViscousForceOnSolid : public FluidViscousForceOnSolid
		{
		public:
			FluidAngularConservativeViscousForceOnSolid(SPHBodyContactRelation* body_contact_relation)
				: FluidViscousForceOnSolid(body_contact_relation) {};
			virtual ~FluidAngularConservativeViscousForceOnSolid() {};
		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FluidPressureForceOnSolid
		* @brief Computing the pressure force from the fluid.
		* The pressrue force is added on the viscous force of the latter is computed.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class FluidPressureForceOnSolid : public FluidViscousForceOnSolid
		{
		public:
			FluidPressureForceOnSolid(SPHBodyContactRelation* body_contact_relation);
			virtual ~FluidPressureForceOnSolid() {};
		protected:
			StdLargeVec<Vecd>& force_from_fluid_, & dvel_dt_ave_;
			StdVec<StdLargeVec<Real>*> contact_p_;
			StdVec<StdLargeVec<Vecd>*> contact_dvel_dt_others_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class TotalViscousForceOnSolid
		* @brief Computing the total viscous force from fluid
		*/
		class TotalViscousForceOnSolid :
			public ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>,
			public SolidDataDelegateSimple
		{
		public:
			explicit TotalViscousForceOnSolid(SolidBody* body);
			virtual ~TotalViscousForceOnSolid() {};
		protected:
			StdLargeVec<Vecd>& viscous_force_from_fluid_;
			Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TotalForceOnSolid
		 * @brief Computing total force from fluid.
		 */
		class TotalForceOnSolid :
			public ParticleDynamicsReduce<Vecd, ReduceSum<Vecd>>,
			public SolidDataDelegateSimple
		{
		public:
			explicit TotalForceOnSolid(SolidBody* body);
			virtual ~TotalForceOnSolid() {};
		protected:
			StdLargeVec<Vecd>& force_from_fluid_;
			Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class CorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class CorrectConfiguration :
			public InteractionDynamics, public SolidDataDelegateInner
		{
		public:
			CorrectConfiguration(SPHBodyInnerRelation* body_inner_relation);
			virtual ~CorrectConfiguration() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Matd>& B_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**@class ConstrainSolidBodyRegion
	 * @brief Constrain a solid body part with prescribed motion.
	 * Note the average values for FSI are prescirbed also.
	 */
		class ConstrainSolidBodyRegion :
			public PartDynamicsByParticle, public SolidDataDelegateSimple
		{
		public:
			ConstrainSolidBodyRegion(SolidBody* body, BodyPartByParticle* body_part);
			virtual ~ConstrainSolidBodyRegion() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, & pos_0_;
			StdLargeVec<Vecd>& n_, & n_0_;
			StdLargeVec<Vecd>& vel_n_, & dvel_dt_, & vel_ave_, & dvel_dt_ave_;
			virtual Point GetDisplacement(Point& pos_0, Point& pos_n) { return pos_0; };
			virtual Vecd GetVelocity(Point& pos_0, Point& pos_n, Vecd& vel_n) { return Vecd(0); };
			virtual Vecd GetAcceleration(Point& pos_0, Point& pos_n, Vecd& dvel_dt) { return Vecd(0); };
			virtual SimTK::Rotation getBodyRotation(Point& pos_0, Point& pos_n, Vecd& dvel_dt) { return SimTK::Rotation(); }
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**@class ImposeExternalForce
		 * @brief impose external force on a solid body part
		 * by add extra acceleration
		 */
		class ImposeExternalForce :
			public PartDynamicsByParticle, public SolidDataDelegateSimple
		{
		public:
			ImposeExternalForce(SolidBody* body, SolidBodyPartForSimbody* body_part);
			virtual ~ImposeExternalForce() {};
		protected:
			StdLargeVec<Vecd>& pos_0_, & vel_n_, & vel_ave_;
			/**
			 * @brief acceleration will be specified by the application
			 */
			virtual Vecd GetAcceleration(Vecd& pos) = 0;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ConstrainSolidBodyPartBySimBody
		 * @brief Constrain a solid body part from the motion
		 * computed from Simbody.
		 */
		class ConstrainSolidBodyPartBySimBody : public ConstrainSolidBodyRegion
		{
		public:
			ConstrainSolidBodyPartBySimBody(SolidBody* body,
				SolidBodyPartForSimbody* body_part,
				SimTK::MultibodySystem& MBsystem,
				SimTK::MobilizedBody& mobod,
				SimTK::Force::DiscreteForces& force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator& integ);
			virtual ~ConstrainSolidBodyPartBySimBody() {};
		protected:
			SimTK::MultibodySystem& MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::Force::DiscreteForces& force_on_bodies_;
			SimTK::RungeKuttaMersonIntegrator& integ_;
			const SimTK::State* simbody_state_;
			Vec3d initial_mobod_origin_location_;

			virtual void setupDynamics(Real dt = 0.0) override;
			void virtual Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TotalForceOnSolidBodyPartForSimBody
		 * @brief Compute the force acting on the solid body part
		 * for applying to simbody forces latter
		 */
		class TotalForceOnSolidBodyPartForSimBody :
			public PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>,
			public SolidDataDelegateSimple
		{
		public:
			TotalForceOnSolidBodyPartForSimBody(SolidBody* body,
				SolidBodyPartForSimbody* body_part,
				SimTK::MultibodySystem& MBsystem,
				SimTK::MobilizedBody& mobod,
				SimTK::Force::DiscreteForces& force_on_bodies,
				SimTK::RungeKuttaMersonIntegrator& integ);
			virtual ~TotalForceOnSolidBodyPartForSimBody() {};
		protected:
			StdLargeVec<Vecd>& force_from_fluid_, & contact_force_, & pos_n_;
			SimTK::MultibodySystem& MBsystem_;
			SimTK::MobilizedBody& mobod_;
			SimTK::Force::DiscreteForces& force_on_bodies_;
			SimTK::RungeKuttaMersonIntegrator& integ_;
			const SimTK::State* simbody_state_;
			Vec3d current_mobod_origin_location_;

			virtual void SetupReduce() override;
			virtual SimTK::SpatialVec ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class FluidViscousForceOnSolid
		* @brief Computing the viscous force from the fluid with wall modeling
		*/
		class FluidViscousForceOnSolidWallModel : public FluidViscousForceOnSolid
		{
		public:
			FluidViscousForceOnSolidWallModel(SPHBodyContactRelation* body_contact_relation,
				fluid_dynamics::ViscousAccelerationWallModel* viscous_acceleration_wall_modeling);
			virtual ~FluidViscousForceOnSolidWallModel() {};
		protected:
			StdLargeVec<Vecd>& n_, & gradient_p_;
			StdLargeVec<Matd>& gradient_vel_;

			//dynamics of a particle
			//to be realized in specific algorithms
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
		
		//----------------------------------------------------------------------
		//		for elastic solid dynamics 
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataDelegateSimple;
		typedef DataDelegateInner<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataDelegateInner;

		/**
		 * @class ElasticSolidDynamicsInitialCondition
		 * @brief  set initial condition for a solid body with different material
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElasticSolidDynamicsInitialCondition : 
			public ParticleDynamicsSimple, public ElasticSolidDataDelegateSimple
		{
		public:
			ElasticSolidDynamicsInitialCondition(SolidBody *body);
			virtual ~ElasticSolidDynamicsInitialCondition() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, & vel_n_;
		};

		/**
		* @class UpdateElasticNormalDirection
		* @brief update particle normal directions for elastic solid
		*/
		class UpdateElasticNormalDirection :
			public ParticleDynamicsSimple, public ElasticSolidDataDelegateSimple
		{
		public:
			explicit UpdateElasticNormalDirection(SolidBody *elastic_body);
			virtual ~UpdateElasticNormalDirection() {};
		protected:
			StdLargeVec<Vecd>& n_, & n_0_;
			StdLargeVec<Matd>& F_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class AcousticTimeStepSize :
			public ParticleDynamicsReduce<Real, ReduceMin>,
			public ElasticSolidDataDelegateSimple
		{
		public:
			explicit AcousticTimeStepSize(SolidBody* body);
			virtual ~AcousticTimeStepSize() {};
		protected:
			StdLargeVec<Vecd>& vel_n_, & dvel_dt_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DeformationGradientTensorBySummation
		* @brief computing deformation gradient tensor by summation
		*/
		class DeformationGradientTensorBySummation :
			public InteractionDynamics, public ElasticSolidDataDelegateInner
		{
		public:
			DeformationGradientTensorBySummation(SPHBodyInnerRelation* body_inner_relation);
			virtual ~DeformationGradientTensorBySummation() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& pos_n_;
			StdLargeVec<Matd>& B_, & F_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class StressRelaxationFirstHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class StressRelaxationFirstHalf 
			: public ParticleDynamics1Level, public ElasticSolidDataDelegateInner
		{
		public:
			StressRelaxationFirstHalf(SPHBodyInnerRelation* body_inner_relation);
			virtual ~StressRelaxationFirstHalf() {};
		protected:
			Real rho_0_, inv_rho_0_;
			StdLargeVec<Real>& Vol_, & rho_n_, & mass_;
			StdLargeVec<Vecd>& pos_n_, & vel_n_, & dvel_dt_, & dvel_dt_others_, & force_from_fluid_;
			StdLargeVec<Matd>& B_, & F_, & dF_dt_, & stress_;
			Real numerical_viscosity_;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class StressRelaxationSecondHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class StressRelaxationSecondHalf : public StressRelaxationFirstHalf
		{
		public:
			StressRelaxationSecondHalf(SPHBodyInnerRelation* body_inner_relation) :
				StressRelaxationFirstHalf(body_inner_relation) {};
			virtual ~StressRelaxationSecondHalf() {};
		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class InitializeDisplacement
		* @brief initialize the displacement for computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class InitializeDisplacement :
			public ParticleDynamicsSimple, public ElasticSolidDataDelegateSimple
		{
		public:
			explicit InitializeDisplacement(SolidBody* body, StdLargeVec<Vecd>& pos_temp);
			virtual ~InitializeDisplacement() {};
		protected:
			StdLargeVec<Vecd>& pos_temp_, & pos_n_, & vel_ave_, & dvel_dt_ave_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class UpdateAverageVelocityAndAcceleration
		* @brief Computing average velocity.
		* This class is for FSI applications to achieve smaller solid dynamics
		* time step size compared to the fluid dynamics
		*/
		class UpdateAverageVelocityAndAcceleration : public InitializeDisplacement
		{
		public:
			explicit UpdateAverageVelocityAndAcceleration(SolidBody* body, StdLargeVec<Vecd>& pos_temp)
				: InitializeDisplacement(body, pos_temp) {};
			virtual ~UpdateAverageVelocityAndAcceleration() {};
		protected:
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AverageVelocityAndAcceleration
		* @brief impose force matching between fluid and solid dynamics
		*/
		class AverageVelocityAndAcceleration
		{
		public:
			InitializeDisplacement initialize_displacement_;
			UpdateAverageVelocityAndAcceleration update_averages_;

			AverageVelocityAndAcceleration(SolidBody* body);
			~AverageVelocityAndAcceleration() {};
		protected:
			StdLargeVec<Vecd> pos_temp_;
		};
	}
}
