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
* @file 	relax_dynamics.h
* @brief 	This is the classes of particle relaxation in order to produce body fitted
* 			initial particle distribution.   
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "math.h"

#include "all_particle_dynamics.h"
#include "neighbor_relation.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"

namespace SPH
{
	namespace relax_dynamics
	{
		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsSimple = ParticleDynamicsSimple<SPHBody, ParticlesType>;

		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsMin = ParticleDynamicsReduce<Real, ReduceMin, SPHBody, ParticlesType>;

		template <class ReturnType, class ParticlesType = BaseParticles>
		using RelaxDynamicsSum = ParticleDynamicsReduce<ReturnType, ReduceSum<ReturnType>, SPHBody,	ParticlesType>;

		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsInner1Level = ParticleDynamicsInner1Level<SPHBody, ParticlesType>;
	
		template <class ParticlesType = BaseParticles>
		using RelaxDynamicsInner =  ParticleDynamicsInner<SPHBody, ParticlesType>;

		template <class ParticlesType = BaseParticles>
		using RelaxSurfaceConstraint = PartDynamicsByParticle<SPHBody, ParticlesType, BodySurface>;

		template <class ParticlesType = BaseParticles>
		using RelaxConstraintByCell = PartDynamicsByCell<SPHBody, ParticlesType, NearBodySurface>;

		template <class ParticlesType = BaseParticles, class ContactParticlesType = BaseParticles>
		using RelaxDynamicsComplex =  
			ParticleDynamicsComplex<SPHBody, ParticlesType, BaseMaterial, SPHBody, ParticlesType>;

		template <class ParticlesType = BaseParticles, class ContactParticlesType = BaseParticles>
		using RelaxDynamicsComplex1Level 
			= ParticleDynamicsComplex1Level<SPHBody, ParticlesType, BaseMaterial, SPHBody, ParticlesType>;
	
		/**
		* @class GetTimeStepSize
		* @brief relaxation dynamics for particle initialization
		* computing time step size
		*/
		class GetTimeStepSize : public RelaxDynamicsMin<BaseParticles>
		{
		protected:
			Real smoothing_length_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit GetTimeStepSize(SPHBody* body);
			virtual ~GetTimeStepSize() {};
		};

		/**
		* @class PhysicsRelaxationInner
		* @brief simple algorithm for physics relaxation
		* without considering contact interaction.
		* this is usually used for solid like bodies
		*/
		class PhysicsRelaxationInner : public ParticleDynamicsInner<SPHBody>
		{
		protected:
			Real smoothing_length_;
			Real particle_spacing_;
			Real sound_speed_;
			Real p0_;
			Real p_star_;
			Real mass_;

			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PhysicsRelaxationInner(SPHBodyInnerRelation* body_inner_relation);
			virtual ~PhysicsRelaxationInner() {};
		};
		/**
		* @class PhysicsRelaxationComplex
		* @brief position verlet algorithm for physics relaxation
		* with considering contact interaction
		* this is usually used for fluid like bodies
		*/
		class PhysicsRelaxationComplex : public RelaxDynamicsComplex1Level<BaseParticles>
		{
		protected:
			Real eta_;
			Real p0_;
			Real sound_speed_;
			Real p_star_;
			Real mass_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			PhysicsRelaxationComplex(SPHBodyComplexRelation* body_complex_relation);
			virtual ~PhysicsRelaxationComplex() {};
		};

		/**
		* @class BodySurfaceBounding
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class BodySurfaceBounding : public RelaxConstraintByCell<BaseParticles>
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			BodySurfaceBounding(SPHBody *body, NearBodySurface* body_part)
				:RelaxConstraintByCell<BaseParticles>(body, body_part) {};
			virtual ~BodySurfaceBounding() {};
		};

		/**
		* @class ConstraintSurfaceParticles
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstraintSurfaceParticles : public RelaxSurfaceConstraint<BaseParticles>
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ConstraintSurfaceParticles(SPHBody* body, BodySurface* body_part)
				:RelaxSurfaceConstraint<BaseParticles>(body, body_part) {};
			virtual ~ConstraintSurfaceParticles() {};
		};

		/**
		* @class computeNumberDensityBySummation
		* @brief  compute the particle number density by summation.
		*/
		class computeNumberDensityBySummation : public RelaxDynamicsComplex<BaseParticles>
		{
		protected:
			Real W0_;
			virtual void ComplexInteraction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			computeNumberDensityBySummation(SPHBodyComplexRelation* body_complex_relation)
				: RelaxDynamicsComplex<BaseParticles>(body_complex_relation) {
				W0_ = body_->kernel_->W(Vecd(0));
			};
			virtual ~computeNumberDensityBySummation() {};
		};

		/**
		 * @class getAveragedParticleNumberDensity
		 * @brief  Compute the Everaged Particle Number Density.
		 */
		class getAveragedParticleNumberDensity  : public RelaxDynamicsSum<Real, BaseParticles>
		{
		protected:
			Real average_farctor_;
			Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit getAveragedParticleNumberDensity(SPHBody* body);
			virtual ~getAveragedParticleNumberDensity() {};
		};
		/**
		* @class FinalizingParticleRelaxation
		* @brief update the number density after relaxation.
		*/
		class FinalizingParticleRelaxation : public RelaxDynamicsSimple<BaseParticles>
		{
		protected:
			/** the average particle number density. */
			Real sigma_;
			/** the method to compute average particle number density. */
			getAveragedParticleNumberDensity* get_average_number_density_;

			/** the function for set global parameters for the particle dynamics */
			virtual void setupDynamics(Real dt = 0.0) override {
				body_->setNewlyUpdated();
				sigma_ = get_average_number_density_->parallel_exec(); 
			};
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			explicit FinalizingParticleRelaxation(SPHBody *body)
				: RelaxDynamicsSimple<BaseParticles>(body), sigma_(0.0){
				get_average_number_density_ = new getAveragedParticleNumberDensity(body);
			};
			virtual ~FinalizingParticleRelaxation() {};
		};
	}
}