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
#include "thin_structure_math.h"

namespace SPH
{
	namespace thin_structure_dynamics
	{
		//----------------------------------------------------------------------
		//		for shell solid dynamics 
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<ThinStructure, ShellParticles, ElasticSolid> ShellDataDelegateSimple;
		typedef DataDelegateInner<ThinStructure, ShellParticles, ElasticSolid> ShellDataDelegateInner;

		/**
		 * @class ShellDynamicsInitialCondition
		 * @brief  set initial condition for shell particles
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ShellDynamicsInitialCondition :
			public ParticleDynamicsSimple, public ShellDataDelegateSimple
		{
		public:
			ShellDynamicsInitialCondition(SolidBody *body);
			virtual ~ShellDynamicsInitialCondition() {};
		protected:
			StdLargeVec<Vecd>& n_0_, &n_, &pseudo_n_, &pos_0_;
		};

		/**
		* @class ShellAcousticTimeStepSize
		* @brief Computing the acoustic time step size for shell
		*/
		class ShellAcousticTimeStepSize :
			public ParticleDynamicsReduce<Real, ReduceMin>,
			public ShellDataDelegateSimple
		{
		public:
			explicit ShellAcousticTimeStepSize(SolidBody* body);
			virtual ~ShellAcousticTimeStepSize() {};
		protected:
			StdLargeVec<Vecd>& vel_n_, &dvel_dt_, &angular_vel_, &dangular_vel_dt_;
			StdLargeVec<Real>& shell_thickness_;
			Real rho_0_, physical_viscosity_, E_0_, nu_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellCorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class ShellCorrectConfiguration :
			public InteractionDynamics, public ShellDataDelegateInner
		{
		public:
			ShellCorrectConfiguration(BaseInnerBodyRelation* body_inner_relation);
			virtual ~ShellCorrectConfiguration() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Matd>& B_;
			StdLargeVec<Vecd>& n_0_;
			StdLargeVec<Matd>& transformation_matrix_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellDeformationGradientTensor
		* @brief computing deformation gradient tensor for shell
		*/
		class ShellDeformationGradientTensor :
			public InteractionDynamics, public ShellDataDelegateInner
		{
		public:
			ShellDeformationGradientTensor(BaseInnerBodyRelation* body_inner_relation);
			virtual ~ShellDeformationGradientTensor() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& pos_n_, &pseudo_n_, &n_0_;
			StdLargeVec<Matd>& B_, &F_, &F_bending_;
			StdLargeVec<Matd>& transformation_matrix_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellStressRelaxationFirstHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class ShellStressRelaxationFirstHalf :
			public ParticleDynamics1Level, public ShellDataDelegateInner
		{
		public:
			ShellStressRelaxationFirstHalf(BaseInnerBodyRelation* body_inner_relation);
			virtual ~ShellStressRelaxationFirstHalf() {};
		protected:
			Real rho_0_, inv_rho_0_;
			StdLargeVec<Real>& Vol_, &rho_n_, &mass_, &shell_thickness_;
			StdLargeVec<Vecd>& pos_n_, &vel_n_, &dvel_dt_, &dvel_dt_others_, &force_from_fluid_;
			StdLargeVec<Vecd>& n_0_, &pseudo_n_, &dpseudo_n_dt_, &dpseudo_n_d2t_, &rotation_, 
				&angular_vel_, dangular_vel_dt_;
			StdLargeVec<Vecd>& shear_stress_;
			StdLargeVec<Matd>& B_, &stress_PK1_, &F_, &dF_dt_, &F_bending_, &dF_bending_dt_,
				&corrected_stress_, &corrected_moment_;
			StdLargeVec<Matd>& transformation_matrix_;
			Real numerical_viscosity_;

			Real shear_correction_factor;
			int number_of_gaussian_point;
			const vector<Real>* gaussian_point;
			const vector<Real>* gaussian_weight;
			static vector<Real> three_gaussian_points;
			static vector<Real> three_gaussian_weights;
			static vector<Real> five_gaussian_points;
			static vector<Real> five_gaussian_weights;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShellStressRelaxationSecondHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class ShellStressRelaxationSecondHalf : public ShellStressRelaxationFirstHalf
		{
		public:
			ShellStressRelaxationSecondHalf(BaseInnerBodyRelation* body_inner_relation)
				: ShellStressRelaxationFirstHalf(body_inner_relation) {};
			virtual ~ShellStressRelaxationSecondHalf() {};
		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**@class ConstrainShellBodyRegion
		 * @brief Fix the position and angle of a shell body part.
		 * Note that the average values for FSI are prescirbed also.
		 */
		class ConstrainShellBodyRegion :
			public PartSimpleDynamicsByParticle, public ShellDataDelegateSimple
		{
		public:
			ConstrainShellBodyRegion(SolidBody* body, BodyPartByParticle* body_part);
			virtual ~ConstrainShellBodyRegion() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, &pos_0_;
			StdLargeVec<Vecd>& n_;
			StdLargeVec<Vecd>& vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			StdLargeVec<Vecd>& rotation_, &angular_vel_, &dangular_vel_dt_;
			StdLargeVec<Vecd>& pseudo_n_, &dpseudo_n_dt_;
			virtual Vecd getDisplacement(Vecd& pos_0, Vecd& pos_n) { return pos_0; };
			virtual Vecd getVelocity(Vecd& pos_0, Vecd& pos_n, Vecd& vel_n) { return Vecd(0); };
			virtual Vecd GetAcceleration(Vecd& pos_0, Vecd& pos_n, Vecd& dvel_dt) { return Vecd(0); };
			virtual Vecd GetRotationAngle(Vecd& pos_0, Vecd& pos_n, Vecd& rotation_angles_0_) { return rotation_angles_0_; };
			virtual Vecd GetAngularVelocity(Vecd& pos_0, Vecd& pos_n, Vecd& angular_vel_) { return Vecd(0); };
			virtual Vecd GetAngularAcceleration(Vecd& pos_0, Vecd& pos_n, Vecd& dangular_vel_dt_) { return Vecd(0); };
			virtual Vecd GetPseudoNormal(Vecd& pos_0, Vecd& pos_n, Vecd& n_0) { return n_0; };
			virtual Vecd GetPseudoNormalChangeRate(Vecd& pos_0, Vecd& pos_n, Vecd& dpseudo_normal_dt_) { return Vecd(0); };
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class FixedFreeRotateShellBoundary
		 * @brief Soft the constrain of a solid body part
		 */
		class FixedFreeRotateShellBoundary :
			public PartInteractionDynamicsByParticle1Level,
			public ShellDataDelegateInner
		{
		public:
			FixedFreeRotateShellBoundary(BaseInnerBodyRelation* body_inner_relation, BodyPartByParticle* body_part);
			virtual ~FixedFreeRotateShellBoundary() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& vel_n_, & angular_vel_;
			StdLargeVec<Vecd>vel_n_temp_, angular_vel_temp_;
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ClampConstrainShellBodyRegion
		 * @brief The clamped constrain of a shell body part
		 */
		class ClampConstrainShellBodyRegion :
			public PartInteractionDynamicsByParticle1Level,
			public ShellDataDelegateInner
		{
		public:
			ClampConstrainShellBodyRegion(BaseInnerBodyRelation* body_inner_relation, BodyPartByParticle* body_part);
			virtual ~ClampConstrainShellBodyRegion() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& vel_n_, &angular_vel_;
			StdLargeVec<Vecd>vel_n_temp_, angular_vel_temp_;
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**@class ConstrainShellBodyRegionInAxisDirection
		 * @brief The boundary conditions are denoted by SS1 according to the references.
	     * The axis_direction must be 0 or 1.
		 * Note that the average values for FSI are prescirbed also.
		 */
		class ConstrainShellBodyRegionInAxisDirection :
			public PartSimpleDynamicsByParticle, public ShellDataDelegateSimple
		{
		public:
			ConstrainShellBodyRegionInAxisDirection(SolidBody* body, BodyPartByParticle* body_part, int axis_direction);
			virtual ~ConstrainShellBodyRegionInAxisDirection() {};
		protected:
			const int axis_;    /**< the axis direction for bounding*/
			StdLargeVec<Vecd>& pos_n_, &pos_0_;
			StdLargeVec<Vecd>& vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			StdLargeVec<Vecd>& rotation_, &angular_vel_, &dangular_vel_dt_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //THIN_STRUCTURE_DYNAMICS_H