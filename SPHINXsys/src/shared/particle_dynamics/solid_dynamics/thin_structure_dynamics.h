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
* @file 	thin_structure_dynamics.h
* @brief 	Here, we define the algorithm classes for thin structure dynamics. 
* @details 	We consider here a weakly compressible solids.   
* @author	Dong Wu, Chi Zhang and Xiangyu Hu
*/

#ifndef THIN_STRUCTURE_DYNAMICS_H
#define THIN_STRUCTURE_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "weakly_compressible_fluid.h"
#include "base_kernel.h"
#include "all_fluid_dynamics.h"

namespace SPH
{
	namespace thin_structure_dynamics
	{
		typedef DataDelegateSimple<ThinStructure, ShellParticles, ElasticSolid> ShellDataSimple;
		typedef DataDelegateInner<ThinStructure, ShellParticles, ElasticSolid> ShellDataInner;

		/**
		 * @class ShellDynamicsInitialCondition
		 * @brief  set initial condition for shell particles
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ShellDynamicsInitialCondition : public ParticleDynamicsSimple, public ShellDataSimple
		{
		public:
			explicit ShellDynamicsInitialCondition(SolidBody &solid_body);
			virtual ~ShellDynamicsInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &n_0_, &n_, &pseudo_n_, &pos_0_;
			StdLargeVec<Matd> &transformation_matrix_;
		};

		/**
		* @class ShellAcousticTimeStepSize
		* @brief Computing the acoustic time step size for shell
		*/
		class ShellAcousticTimeStepSize : public ParticleDynamicsReduce<Real, ReduceMin>,
										  public ShellDataSimple
		{
		public:
			explicit ShellAcousticTimeStepSize(SolidBody &sph_body);
			virtual ~ShellAcousticTimeStepSize(){};

		protected:
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &angular_vel_, &dangular_vel_dt_;
			StdLargeVec<Real> &shell_thickness_;
			Real rho0_, physical_viscosity_, E0_, nu_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellCorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class ShellCorrectConfiguration : public InteractionDynamics, public ShellDataInner
		{
		public:
			explicit ShellCorrectConfiguration(BaseBodyRelationInner &inner_relation);
			virtual ~ShellCorrectConfiguration(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Matd> &B_;
			StdLargeVec<Vecd> &n_0_;
			StdLargeVec<Matd> &transformation_matrix_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellDeformationGradientTensor
		* @brief computing deformation gradient tensor for shell
		*/
		class ShellDeformationGradientTensor : public InteractionDynamics, public ShellDataInner
		{
		public:
			explicit ShellDeformationGradientTensor(BaseBodyRelationInner &inner_relation);
			virtual ~ShellDeformationGradientTensor(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &pos_n_, &pseudo_n_, &n_0_;
			StdLargeVec<Matd> &B_, &F_, &F_bending_;
			StdLargeVec<Matd> &transformation_matrix_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class BaseShellRelaxation
		 * @brief abstract class for preparing shell relaxation
		*/
		class BaseShellRelaxation : public ParticleDynamics1Level, public ShellDataInner
		{
		public:
			explicit BaseShellRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseShellRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &rho_n_, &mass_, &shell_thickness_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dvel_dt_, &dvel_dt_prior_, &force_from_fluid_;
			StdLargeVec<Vecd> &n_0_, &pseudo_n_, &dpseudo_n_dt_, &dpseudo_n_d2t_, &rotation_,
				&angular_vel_, dangular_vel_dt_;
			StdLargeVec<Matd> &B_, &F_, &dF_dt_, &F_bending_, &dF_bending_dt_;
			StdLargeVec<Matd> &transformation_matrix_;
		};

		/**
		* @class ShellStressRelaxationFirstHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class ShellStressRelaxationFirstHalf : public BaseShellRelaxation
		{
		public:
			explicit ShellStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation,
				int number_of_gaussian_points = 3, bool hourglass_control = false);
			virtual ~ShellStressRelaxationFirstHalf() {};

		protected:
			Real rho0_, inv_rho0_;
			StdLargeVec<Matd> &stress_PK1_, &global_stress_, &global_moment_;
			StdLargeVec<Vecd> &global_shear_stress_, &n_;
			Real smoothing_length_, E0_, G0_, nu_, hourglass_control_factor_;
			bool hourglass_control_;
			const Real inv_W0_ = 1.0 / body_->sph_adaptation_->getKernel()->W0(Vecd(0));
			const Real shear_correction_factor_ = 5.0 / 6.0;

			const StdVec<Real> three_gaussian_points_ = { 0.0, 0.7745966692414834, -0.7745966692414834 };
			const StdVec<Real> three_gaussian_weights_ = { 0.8888888888888889, 0.5555555555555556, 0.5555555555555556 };
			const StdVec<Real> five_gaussian_points_
				= { 0.0, 0.5384693101056831, -0.5384693101056831, 0.9061798459386640, -0.9061798459386640 };
			const StdVec<Real> five_gaussian_weights_
				= { 0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891 };
			int number_of_gaussian_points_;
			StdVec<Real> gaussian_point_;
			StdVec<Real> gaussian_weight_;


			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellStressRelaxationSecondHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class ShellStressRelaxationSecondHalf : public BaseShellRelaxation
		{
		public:
			explicit ShellStressRelaxationSecondHalf(BaseBodyRelationInner &inner_relation)
				: BaseShellRelaxation(inner_relation){};
			virtual ~ShellStressRelaxationSecondHalf(){};

		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**@class ConstrainShellBodyRegion
		 * @brief Fix the position and angle of a shell body part.
		 * Note that the average values for FSI are prescirbed also.
		 */
		class ConstrainShellBodyRegion : public PartSimpleDynamicsByParticle, public ShellDataSimple
		{
		public:
			ConstrainShellBodyRegion(SolidBody &sph_body, BodyPartByParticle &body_part);
			virtual ~ConstrainShellBodyRegion(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> &n_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			StdLargeVec<Vecd> &rotation_, &angular_vel_, &dangular_vel_dt_;
			StdLargeVec<Vecd> &pseudo_n_, &dpseudo_n_dt_;
			virtual Vecd getDisplacement(const Vecd &pos_0, const Vecd &pos_n) { return pos_0; };
			virtual Vecd getVelocity(const Vecd &pos_0, const Vecd &pos_n, const Vecd &vel_n) { return Vecd(0); };
			virtual Vecd GetAcceleration(const Vecd &pos_0, const Vecd &pos_n, const Vecd &dvel_dt) { return Vecd(0); };
			virtual Vecd GetRotationAngle(const Vecd &pos_0, const Vecd &pos_n, const Vecd &rotation_angles_0_) { return rotation_angles_0_; };
			virtual Vecd GetAngularVelocity(const Vecd &pos_0, const Vecd &pos_n, const Vecd &angular_vel_) { return Vecd(0); };
			virtual Vecd GetAngularAcceleration(const Vecd &pos_0, const Vecd &pos_n, const Vecd &dangular_vel_dt_) { return Vecd(0); };
			virtual Vecd GetPseudoNormal(const Vecd &pos_0, const Vecd &pos_n, const Vecd &n_0) { return n_0; };
			virtual Vecd GetPseudoNormalChangeRate(const Vecd &pos_0, const Vecd &pos_n, const Vecd &dpseudo_normal_dt_) { return Vecd(0); };
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class FixedFreeRotateShellBoundary
		 * @brief Soft the constraint of a solid body part
		 */
		class FixedFreeRotateShellBoundary : public PartInteractionDynamicsByParticle1Level,
											 public ShellDataInner
		{
		public:
			FixedFreeRotateShellBoundary(BaseBodyRelationInner &inner_relation,
										 BodyPartByParticle &body_part, Vecd constrained_direction = Vecd(0));
			virtual ~FixedFreeRotateShellBoundary(){};

		protected:
			Real W0_;
			Matd constrain_matrix_, recover_matrix_;
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_n_, &angular_vel_;
			StdLargeVec<Vecd> vel_n_temp_, angular_vel_temp_;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ClampConstrainShellBodyRegion
		 * @brief The clamped constrain of a shell body part
		 */
		class ClampConstrainShellBodyRegion : public PartInteractionDynamicsByParticle1Level,
											  public ShellDataInner
		{
		public:
			ClampConstrainShellBodyRegion(BaseBodyRelationInner &inner_relation, BodyPartByParticle &body_part);
			virtual ~ClampConstrainShellBodyRegion(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_n_, &angular_vel_;
			StdLargeVec<Vecd> vel_n_temp_, angular_vel_temp_;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**@class ConstrainShellBodyRegionInAxisDirection
		 * @brief The boundary conditions are denoted by SS1 according to the references.
	     * The axis_direction must be 0 or 1.
		 * Note that the average values for FSI are prescirbed also.
		 */
		class ConstrainShellBodyRegionInAxisDirection : public PartSimpleDynamicsByParticle, public ShellDataSimple
		{
		public:
			ConstrainShellBodyRegionInAxisDirection(SolidBody &sph_body, BodyPartByParticle &body_part, int axis_direction);
			virtual ~ConstrainShellBodyRegionInAxisDirection(){};

		protected:
			const int axis_; /**< the axis direction for bounding*/
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			StdLargeVec<Vecd> &rotation_, &angular_vel_, &dangular_vel_dt_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class DistributingPointForcesToShell
		 * @brief Distribute a series of point forces to its contact shell bodies.
		 */
		class DistributingPointForcesToShell : public ParticleDynamicsSimple, public ShellDataSimple
		{
		protected:
			std::vector<Vecd> point_forces_, reference_positions_, time_dependent_point_forces_;
			Real time_to_full_external_force_;
			Real particle_spacing_ref_, h_spacing_ratio_;
			StdLargeVec<Vecd> &pos_0_, &dvel_dt_prior_;
			StdLargeVec<Real> &Vol_, &mass_, &shell_thickness_;
			std::vector <StdLargeVec<Real>> weight_;
			std::vector<Real> sum_of_weight_;

		public:
			DistributingPointForcesToShell(SolidBody &sph_body, std::vector<Vecd> point_forces,
				std::vector<Vecd> reference_positions, Real time_to_full_external_force,
				Real particle_spacing_ref, Real h_spacing_ratio = 1.6);
			virtual ~DistributingPointForcesToShell() {};

			void getWeight();
			void getForce();
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //THIN_STRUCTURE_DYNAMICS_H