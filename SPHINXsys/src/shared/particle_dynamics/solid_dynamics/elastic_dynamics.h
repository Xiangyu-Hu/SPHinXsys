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
* @file 	elastic_dynamics.h
* @brief 	Here, we define the algorithm classes for elastic solid dynamics. 
* @details 	We consider here a weakly compressible solids.   
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
*/

#ifndef ELASTIC_DYNAMICS_H
#define ELASTIC_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "base_kernel.h"
#include "body_relation.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"



namespace SPH
{
	template <typename VariableType>
	class BodySummation;
	template <typename VariableType>
	class BodyMoment;

	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for elastic solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataSimple;
		typedef DataDelegateInner<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataInner;

		/**
		 * @class ElasticDynamicsInitialCondition
		 * @brief  set initial condition for a solid body with different material
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElasticDynamicsInitialCondition : public ParticleDynamicsSimple, public ElasticSolidDataSimple
		{
		public:
			explicit ElasticDynamicsInitialCondition(SolidBody &solid_body);
			virtual ~ElasticDynamicsInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_;
		};

		/**
		* @class UpdateElasticNormalDirection
		* @brief update particle normal directions for elastic solid
		*/
		class UpdateElasticNormalDirection : public ParticleDynamicsSimple, public ElasticSolidDataSimple
		{
		public:
			explicit UpdateElasticNormalDirection(SolidBody &solid_body);
			virtual ~UpdateElasticNormalDirection(){};

		protected:
			StdLargeVec<Vecd> &n_, &n_0_;
			StdLargeVec<Matd> &F_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class AcousticTimeStepSize : public ParticleDynamicsReduce<Real, ReduceMin>,
									 public ElasticSolidDataSimple
		{
		public:
			explicit AcousticTimeStepSize(SolidBody &solid_body, Real CFL = 0.6);
			virtual ~AcousticTimeStepSize(){};

		protected:
			Real CFL_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_;
			Real smoothing_length_, c0_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DeformationGradientTensorBySummation
		* @brief computing deformation gradient tensor by summation
		*/
		class DeformationGradientTensorBySummation : public InteractionDynamics, public ElasticSolidDataInner
		{
		public:
			explicit DeformationGradientTensorBySummation(BaseBodyRelationInner &inner_relation);
			virtual ~DeformationGradientTensorBySummation(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &pos_n_;
			StdLargeVec<Matd> &B_, &F_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class BaseElasticRelaxation
		* @brief base class for elastic relaxation
		*/
		class BaseElasticRelaxation
			: public ParticleDynamics1Level,
			  public ElasticSolidDataInner
		{
		public:
			explicit BaseElasticRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseElasticRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &rho_n_, &mass_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dvel_dt_;
			StdLargeVec<Matd> &B_, &F_, &dF_dt_;
		};

		/**
		* @class StressRelaxationFirstHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class StressRelaxationFirstHalf : public BaseElasticRelaxation
		{
		public:
			explicit StressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~StressRelaxationFirstHalf(){};

		protected:
			Real rho0_, inv_rho0_;
			StdLargeVec<Vecd> &dvel_dt_prior_, &force_from_fluid_;
			StdLargeVec<Matd> &stress_PK1_;
			Real numerical_dissipation_factor_;
			Real smoothing_length_;
			Real inv_W0_ = 1.0 / body_->sph_adaptation_->getKernel()->W0(Vecd(0));

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class KirchhoffParticleStressRelaxationFirstHalf
		*/
		class KirchhoffParticleStressRelaxationFirstHalf : public StressRelaxationFirstHalf
		{
		public:
			explicit KirchhoffParticleStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~KirchhoffParticleStressRelaxationFirstHalf(){};

		protected:
			const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class KirchhoffStressRelaxationFirstHalf
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
		class KirchhoffStressRelaxationFirstHalf : public StressRelaxationFirstHalf
		{
		public:
			explicit KirchhoffStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~KirchhoffStressRelaxationFirstHalf(){};

		protected:
			StdLargeVec<Real> J_to_minus_2_over_dimension_;
			StdLargeVec<Matd> stress_on_particle_, inverse_F_T_;
			const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;
			const Real correction_factor_ = 1.07;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class StressRelaxationSecondHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class StressRelaxationSecondHalf : public BaseElasticRelaxation
		{
		public:
			explicit StressRelaxationSecondHalf(BaseBodyRelationInner &inner_relation)
				: BaseElasticRelaxation(inner_relation){};
			virtual ~StressRelaxationSecondHalf(){};

		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
    }
}
#endif //ELASTIC_DYNAMICS_H