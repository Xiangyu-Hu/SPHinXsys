/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	elastic_dynamics.h
 * @brief 	Here, we define the algorithm classes for elastic solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef ELASTIC_DYNAMICS_H
#define ELASTIC_DYNAMICS_H

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
		//		for elastic solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<ElasticSolidParticles> ElasticSolidDataSimple;
		typedef DataDelegateInner<ElasticSolidParticles> ElasticSolidDataInner;

		/**
		 * @class ElasticDynamicsInitialCondition
		 * @brief  set initial condition for a solid body with different material
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElasticDynamicsInitialCondition : public LocalDynamics, public ElasticSolidDataSimple
		{
		public:
			explicit ElasticDynamicsInitialCondition(SPHBody &sph_body);
			virtual ~ElasticDynamicsInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_, &vel_;
		};

		/**
		 * @class UpdateElasticNormalDirection
		 * @brief update particle normal directions for elastic solid
		 */
		class UpdateElasticNormalDirection : public LocalDynamics, public ElasticSolidDataSimple
		{
		protected:
			StdLargeVec<Vecd> &n_, &n0_;
			StdLargeVec<Matd> &F_;

		public:
			explicit UpdateElasticNormalDirection(SPHBody &sph_body);
			virtual ~UpdateElasticNormalDirection(){};

			void update(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class AcousticTimeStepSize
		 * @brief Computing the acoustic time step size
		 * computing time step size
		 */
		class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
									 public ElasticSolidDataSimple
		{
		protected:
			Real CFL_;
			StdLargeVec<Vecd> &vel_, &acc_;
			Real smoothing_length_, c0_;

		public:
			explicit AcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
			virtual ~AcousticTimeStepSize(){};

			Real reduce(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class DeformationGradientBySummation
		 * @brief computing deformation gradient tensor by summation
		 */
		class DeformationGradientBySummation : public LocalDynamics, public ElasticSolidDataInner
		{
		public:
			explicit DeformationGradientBySummation(BaseInnerRelation &inner_relation);
			virtual ~DeformationGradientBySummation(){};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd> &pos_;
			StdLargeVec<Matd> &B_, &F_;
		};

		/**
		 * @class BaseElasticIntegration
		 * @brief base class for elastic relaxation
		 */
		class BaseElasticIntegration : public LocalDynamics, public ElasticSolidDataInner
		{
		public:
			explicit BaseElasticIntegration(BaseInnerRelation &inner_relation);
			virtual ~BaseElasticIntegration(){};

		protected:
			StdLargeVec<Real> &rho_, &mass_;
			StdLargeVec<Vecd> &pos_, &vel_, &acc_;
			StdLargeVec<Matd> &B_, &F_, &dF_dt_;
		};

		/**
		 * @class BaseIntegration1stHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class BaseIntegration1stHalf : public BaseElasticIntegration
		{
		public:
			explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~BaseIntegration1stHalf(){};
			void update(size_t index_i, Real dt = 0.0);

		protected:
			ElasticSolid &elastic_solid_;
			Real rho0_, inv_rho0_;
			StdLargeVec<Vecd> &acc_prior_;
			Real smoothing_length_;
		};

		/**
		 * @class Integration1stHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class Integration1stHalf : public BaseIntegration1stHalf
		{
		public:
			explicit Integration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~Integration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Matd> stress_PK1_B_;
			Real numerical_dissipation_factor_;
			Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(zero_vec);
		};

		/**
		 * @class KirchhoffParticleIntegration1stHalf
		 */
		class KirchhoffParticleIntegration1stHalf : public Integration1stHalf
		{
		public:
			explicit KirchhoffParticleIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~KirchhoffParticleIntegration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
		};

		/**
		 * @class KirchhoffIntegration1stHalf
		 * @brief Decompose the stress into particle stress includes isotropic stress
		 * and the stress due to non-homogeneous material properties.
		 * The preliminary shear stress is introduced by particle pair to avoid
		 * spurious stress and deformation.
		 * Note that, for the shear stress term,
		 * due to the mismatch of the divergence contribution between
		 * the pair-wise second-order derivative Laplacian formulation
		 * and particle-wise first-order gradient formulation,
		 * a correction factor slight large than one is introduced.
		 * Note that, if you see time step size goes unusually small,
		 * it may be due to the determinate of deformation matrix become negative.
		 * In this case, you may need decrease CFL number when computing time-step size.
		 */
		class KirchhoffIntegration1stHalf : public BaseIntegration1stHalf
		{
		public:
			explicit KirchhoffIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~KirchhoffIntegration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Real> J_to_minus_2_over_dimension_;
			StdLargeVec<Matd> stress_on_particle_, inverse_F_T_;
			const Real correction_factor_ = 1.07;
		};

		/**
		 * @class Integration2ndHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the second step
		 */
		class Integration2ndHalf : public BaseElasticIntegration
		{
		public:
			explicit Integration2ndHalf(BaseInnerRelation &inner_relation)
				: BaseElasticIntegration(inner_relation){};
			virtual ~Integration2ndHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif // ELASTIC_DYNAMICS_H