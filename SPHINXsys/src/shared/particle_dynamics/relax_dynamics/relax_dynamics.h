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
	class ComplexShape;

	namespace relax_dynamics
	{
		typedef DataDelegateSimple<SPHBody, BaseParticles> RelaxDataDelegateSimple;

		typedef DataDelegateInner<SPHBody, BaseParticles> RelaxDataDelegateInner;

		typedef DataDelegateComplex<SPHBody, BaseParticles, BaseMaterial, SPHBody, BaseParticles> RelaxDataDelegateComplex;

		/**
		* @class GetTimeStepSizeSquare
		* @brief relaxation dynamics for particle initialization
		* computing the square of time step size
		*/
		class GetTimeStepSizeSquare :
			public ParticleDynamicsReduce<Real, ReduceMin>,
			public RelaxDataDelegateSimple
		{
		public:
			explicit GetTimeStepSizeSquare(SPHBody* body);
			virtual ~GetTimeStepSizeSquare() {};
		protected:
			StdLargeVec<Vecd>& dvel_dt_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class RelaxationAccelerationInner
		* @brief simple algorithm for physics relaxation
		* without considering contact interaction.
		* this is usually used for solid like bodies
		*/
		class RelaxationAccelerationInner : 
			public ParticleDynamicsInner, public RelaxDataDelegateInner
		{
		public:
			RelaxationAccelerationInner(SPHBodyInnerRelation* body_inner_relation);
			virtual ~RelaxationAccelerationInner() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& dvel_dt_;
			StdLargeVec<Vecd>& pos_n_;
			ComplexShape* complex_shape_;
			Kernel* kernel_;
			virtual void InnerInteraction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class UpdateParticlePosition
		* @brief update the particle position for a time step
		*/
		class UpdateParticlePosition :
			public ParticleDynamicsSimple, public RelaxDataDelegateSimple
		{
		public:
			explicit UpdateParticlePosition(SPHBody* body);
			virtual ~UpdateParticlePosition() {};
		protected:
			StdLargeVec<Vecd>& pos_n_, & dvel_dt_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class RelaxationAccelerationComplex
		* @brief compute relaxation acceleration while consider the present of contact bodies
		* with considering contact interaction
		* this is usually used for fluid like bodies
		*/
		class RelaxationAccelerationComplex : 
			public ParticleDynamicsComplex,
			public RelaxDataDelegateComplex
		{
		public:
			RelaxationAccelerationComplex(SPHBodyComplexRelation* body_complex_relation);
			virtual ~RelaxationAccelerationComplex() {};
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Vecd>& dvel_dt_;
			StdVec<StdLargeVec<Real>*> contact_Vol_;
			virtual void ComplexInteraction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class BodySurfaceBounding
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class BodySurfaceBounding : 
			public PartDynamicsByCell,
			public RelaxDataDelegateSimple
		{
		public:
			BodySurfaceBounding(SPHBody *body, NearBodySurface* body_part);
			virtual ~BodySurfaceBounding() {};
		protected:
			StdLargeVec<Vecd>& pos_n_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ConstraintSurfaceParticles
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstraintSurfaceParticles : 
			public PartDynamicsByParticle,
			public RelaxDataDelegateSimple
		{
		public:
			ConstraintSurfaceParticles(SPHBody* body, BodySurface* body_part);
			virtual ~ConstraintSurfaceParticles() {};
		protected:
			StdLargeVec<Vecd>& pos_n_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class RelaxationStepInner
		* @brief carry out particle relaxation step of particles within the body
		*/
		class RelaxationStepInner : public  ParticleDynamics<void>
		{
		protected:
			SPHBody* sph_body_;
			SPHBodyInnerRelation* inner_relation_;
		public:
			explicit RelaxationStepInner(SPHBodyInnerRelation* body_inner_relation);
			virtual ~RelaxationStepInner() {};

			RelaxationAccelerationInner relaxation_acceleration_inner_;
			GetTimeStepSizeSquare get_time_step_square_;
			UpdateParticlePosition update_particle_position_;
			BodySurfaceBounding	surface_bounding_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		* @class computeNumberDensityBySummation
		* @brief  compute the particle number density by summation.
		*/
		class computeNumberDensityBySummation :
			public ParticleDynamicsComplex,
			public RelaxDataDelegateComplex
		{
		public:
			computeNumberDensityBySummation(SPHBodyComplexRelation* body_complex_relation);
			virtual ~computeNumberDensityBySummation() {};
		protected:
			StdLargeVec<Real>& Vol_0_, & sigma_0_;
			StdVec<StdLargeVec<Real>*> contact_Vol_0_;
			Real W0_;
			virtual void ComplexInteraction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class getAveragedParticleNumberDensity
		 * @brief  Compute the Everaged Particle Number Density.
		 */
		class getAveragedParticleNumberDensity : 
			public ParticleDynamicsReduce <Real, ReduceSum<Real>>,
			public RelaxDataDelegateSimple
		{
		public:
			explicit getAveragedParticleNumberDensity(SPHBody* body);
			virtual ~getAveragedParticleNumberDensity() {};
		protected:
			StdLargeVec<Real>& sigma_0_;
			Real average_farctor_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};
		/**
		* @class FinalizingParticleRelaxation
		* @brief update the number density after relaxation.
		*/
		class FinalizingParticleRelaxation : 
			public ParticleDynamicsSimple,
			public RelaxDataDelegateSimple
		{
		public:
			explicit FinalizingParticleRelaxation(SPHBody *body);
			virtual ~FinalizingParticleRelaxation() {};
		protected:
			StdLargeVec<Real>& sigma_0_;
			StdLargeVec<Vecd>& pos_n_;
			/** the average particle number density. */
			Real sigma_;
			/** the method to compute average particle number density. */
			getAveragedParticleNumberDensity* get_average_number_density_;

			/** the function for set global parameters for the particle dynamics */
			virtual void setupDynamics(Real dt = 0.0) override {
				body_->setNewlyUpdated();
				sigma_ = get_average_number_density_->parallel_exec();
			};
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
	}
}