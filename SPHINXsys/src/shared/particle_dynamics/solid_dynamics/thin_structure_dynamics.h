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
 * @file 	thin_structure_dynamics.h
 * @brief 	Here, we define the algorithm classes for thin structure dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Dong Wu, Chi Zhang and Xiangyu Hu
 */

#ifndef THIN_STRUCTURE_DYNAMICS_H
#define THIN_STRUCTURE_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "all_body_relations.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"

namespace SPH
{
	namespace thin_structure_dynamics
	{
		typedef DataDelegateSimple<ShellParticles> ShellDataSimple;
		typedef DataDelegateInner<ShellParticles> ShellDataInner;

		/**
		 * @class ShellDynamicsInitialCondition
		 * @brief  set initial condition for shell particles
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ShellDynamicsInitialCondition : public LocalDynamics, public ShellDataSimple
		{
		public:
			explicit ShellDynamicsInitialCondition(SPHBody &sph_body);
			virtual ~ShellDynamicsInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &n0_, &n_, &pseudo_n_, &pos0_;
			StdLargeVec<Matd> &transformation_matrix_;
		};

		/**
		 * @class ShellAcousticTimeStepSize
		 * @brief Computing the acoustic time step size for shell
		 */
		class ShellAcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
										  public ShellDataSimple
		{
		protected:
			StdLargeVec<Vecd> &vel_, &acc_, &angular_vel_, &dangular_vel_dt_;
			StdLargeVec<Real> &thickness_;
			Real rho0_, physical_viscosity_, E0_, nu_, c0_;
			Real smoothing_length_;
			Real CFL_;

		public:
			explicit ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
			virtual ~ShellAcousticTimeStepSize(){};

			Real reduce(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class ShellCorrectConfiguration
		 * @brief obtain the corrected initial configuration in strong form
		 */
		class ShellCorrectConfiguration : public LocalDynamics, public ShellDataInner
		{
		public:
			explicit ShellCorrectConfiguration(BaseInnerRelation &inner_relation);
			virtual ~ShellCorrectConfiguration(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Matd> &B_;
			StdLargeVec<Vecd> &n0_;
			StdLargeVec<Matd> &transformation_matrix_;
		};

		/**
		 * @class ShellDeformationGradientTensor
		 * @brief computing deformation gradient tensor for shell
		 * TODO: need a test case for this.
		 */
		class ShellDeformationGradientTensor : public LocalDynamics, public ShellDataInner
		{
		public:
			explicit ShellDeformationGradientTensor(BaseInnerRelation &inner_relation);
			virtual ~ShellDeformationGradientTensor(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd> &pos_, &pseudo_n_, &n0_;
			StdLargeVec<Matd> &B_, &F_, &F_bending_;
			StdLargeVec<Matd> &transformation_matrix_;
		};

		/**
		 * @class BaseShellRelaxation
		 * @brief abstract class for preparing shell relaxation
		 */
		class BaseShellRelaxation : public LocalDynamics, public ShellDataInner
		{
		public:
			explicit BaseShellRelaxation(BaseInnerRelation &inner_relation);
			virtual ~BaseShellRelaxation(){};

		protected:
			StdLargeVec<Real> &rho_, &thickness_;
			StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
			StdLargeVec<Vecd> &n0_, &pseudo_n_, &dpseudo_n_dt_, &dpseudo_n_d2t_, &rotation_,
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
			explicit ShellStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
													int number_of_gaussian_points = 3, bool hourglass_control = false);
			virtual ~ShellStressRelaxationFirstHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid &elastic_solid_;
			Real rho0_, inv_rho0_;
			StdLargeVec<Matd> &global_stress_, &global_moment_;
			StdLargeVec<Vecd> &global_shear_stress_, &n_;
			Real smoothing_length_, E0_, G0_, nu_, hourglass_control_factor_;
			bool hourglass_control_;
			const Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(Vecd(0));
			const Real shear_correction_factor_ = 5.0 / 6.0;

			const StdVec<Real> three_gaussian_points_ = {0.0, 0.7745966692414834, -0.7745966692414834};
			const StdVec<Real> three_gaussian_weights_ = {0.8888888888888889, 0.5555555555555556, 0.5555555555555556};
			const StdVec<Real> five_gaussian_points_ = {0.0, 0.5384693101056831, -0.5384693101056831, 0.9061798459386640, -0.9061798459386640};
			const StdVec<Real> five_gaussian_weights_ = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891};
			int number_of_gaussian_points_;
			StdVec<Real> gaussian_point_;
			StdVec<Real> gaussian_weight_;
		};

		/**
		 * @class ShellStressRelaxationSecondHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the second step
		 */
		class ShellStressRelaxationSecondHalf : public BaseShellRelaxation
		{
		public:
			explicit ShellStressRelaxationSecondHalf(BaseInnerRelation &inner_relation)
				: BaseShellRelaxation(inner_relation){};
			virtual ~ShellStressRelaxationSecondHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		};

		/**@class ConstrainShellBodyRegion
		 * @brief Fix the position and angle of a shell body part.
		 * Note that the average values for FSI are prescribed also.
		 */
		class ConstrainShellBodyRegion : public LocalDynamics, public ShellDataSimple
		{
		public:
			ConstrainShellBodyRegion(BodyPartByParticle &body_part);
			virtual ~ConstrainShellBodyRegion(){};
			void update(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd> &pos_, &pos0_;
			StdLargeVec<Vecd> &n_;
			StdLargeVec<Vecd> &vel_, &acc_;
			StdLargeVec<Vecd> &rotation_, &angular_vel_, &dangular_vel_dt_;
			StdLargeVec<Vecd> &pseudo_n_, &dpseudo_n_dt_;
			virtual Vecd getDisplacement(const Vecd &pos_0, const Vecd &pos_n) { return pos_0; };
			virtual Vecd getVelocity(const Vecd &pos_0, const Vecd &pos_n, const Vecd &vel_n) { return Vecd(0); };
			virtual Vecd GetAcceleration(const Vecd &pos_0, const Vecd &pos_n, const Vecd &acc) { return Vecd(0); };
			virtual Vecd GetRotationAngle(const Vecd &pos_0, const Vecd &pos_n, const Vecd &rotation_angles_0_) { return rotation_angles_0_; };
			virtual Vecd GetAngularVelocity(const Vecd &pos_0, const Vecd &pos_n, const Vecd &angular_vel_) { return Vecd(0); };
			virtual Vecd GetAngularAcceleration(const Vecd &pos_0, const Vecd &pos_n, const Vecd &dangular_vel_dt_) { return Vecd(0); };
			virtual Vecd GetPseudoNormal(const Vecd &pos_0, const Vecd &pos_n, const Vecd &n_0) { return n_0; };
			virtual Vecd GetPseudoNormalChangeRate(const Vecd &pos_0, const Vecd &pos_n, const Vecd &dpseudo_normal_dt_) { return Vecd(0); };
		};

		/**@class ConstrainShellBodyRegionAlongAxis
		 * @brief The boundary conditions are denoted by SS1 according to the references.
		 * The axis must be 0 or 1.
		 * Note that the average values for FSI are prescribed also.
		 */
		class ConstrainShellBodyRegionAlongAxis : public LocalDynamics, public ShellDataSimple
		{
		public:
			ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis);
			virtual ~ConstrainShellBodyRegionAlongAxis(){};
			void update(size_t index_i, Real dt = 0.0);

		protected:
			const int axis_; /**< the axis direction for bounding*/
			StdLargeVec<Vecd> &pos_, &pos0_;
			StdLargeVec<Vecd> &vel_, &acc_;
			StdLargeVec<Vecd> &rotation_, &angular_vel_, &dangular_vel_dt_;
		};

		/**
		 * @class DistributingPointForcesToShell
		 * @brief Distribute a series of point forces to its contact shell bodies.
		 */
		class DistributingPointForcesToShell : public LocalDynamics, public ShellDataSimple
		{
		protected:
			std::vector<Vecd> point_forces_, reference_positions_, time_dependent_point_forces_;
			Real time_to_full_external_force_;
			Real particle_spacing_ref_, h_spacing_ratio_;
			StdLargeVec<Vecd> &pos0_, &acc_prior_;
			StdLargeVec<Real> &thickness_;
			std::vector<StdLargeVec<Real>> weight_;
			std::vector<Real> sum_of_weight_;

			void getWeight();

		public:
			DistributingPointForcesToShell(SPHBody &sph_body, std::vector<Vecd> point_forces,
										   std::vector<Vecd> reference_positions, Real time_to_full_external_force,
										   Real particle_spacing_ref, Real h_spacing_ratio = 1.6);
			virtual ~DistributingPointForcesToShell(){};

			virtual void setupDynamics(Real dt = 0.0) override;
			void update(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif // THIN_STRUCTURE_DYNAMICS_H