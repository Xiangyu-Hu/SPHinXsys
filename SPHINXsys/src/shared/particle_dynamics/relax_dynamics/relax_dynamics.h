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
*/

#ifndef RELAX_DYNAMICS_H
#define RELAX_DYNAMICS_H

#include "math.h"

#include "all_particle_dynamics.h"
#include "neighbor_relation.h"
#include "base_kernel.h"
#include "cell_linked_list.h"
#include "solid_dynamics.h"

namespace SPH
{
	class GeometryShape;
	class LevelSetShape;

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
		class GetTimeStepSizeSquare : public ParticleDynamicsReduce<Real, ReduceMax>,
									  public RelaxDataDelegateSimple
		{
		public:
			explicit GetTimeStepSizeSquare(SPHBody &sph_body);
			virtual ~GetTimeStepSizeSquare(){};

		protected:
			StdLargeVec<Vecd> &dvel_dt_;
			Real h_ref_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
			Real OutputResult(Real reduced_value) override;
		};

		/**
		* @class RelaxationAccelerationInner
		* @brief simple algorithm for physics relaxation
		* without considering contact interaction.
		* this is usually used for solid like bodies
		*/
		class RelaxationAccelerationInner : public InteractionDynamics, public RelaxDataDelegateInner
		{
		public:
			explicit RelaxationAccelerationInner(BaseBodyRelationInner &inner_relation);
			virtual ~RelaxationAccelerationInner(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &dvel_dt_;
			StdLargeVec<Vecd> &pos_n_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class RelaxationAccelerationInnerWithLevelSetCorrection
		* @brief we constrain particles to a level function representing the interafce.
		*/
		class RelaxationAccelerationInnerWithLevelSetCorrection : public RelaxationAccelerationInner
		{
		public:
			explicit RelaxationAccelerationInnerWithLevelSetCorrection(BaseBodyRelationInner &inner_relation);
			virtual ~RelaxationAccelerationInnerWithLevelSetCorrection(){};

		protected:
			LevelSetShape *level_set_shape_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class UpdateParticlePosition
		* @brief update the particle position for a time step
		*/
		class UpdateParticlePosition : public ParticleDynamicsSimple,
									   public RelaxDataDelegateSimple
		{
		public:
			explicit UpdateParticlePosition(SPHBody &sph_body);
			virtual ~UpdateParticlePosition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &dvel_dt_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class UpdateSmoothingLengthRatioByBodyShape
		 * @brief update the particle smoothing length ratio
		*/
		class UpdateSmoothingLengthRatioByBodyShape : public ParticleDynamicsSimple,
													  public RelaxDataDelegateSimple
		{
		public:
			explicit UpdateSmoothingLengthRatioByBodyShape(SPHBody &sph_body);
			virtual ~UpdateSmoothingLengthRatioByBodyShape(){};

		protected:
			StdLargeVec<Real> &h_ratio_, &Vol_;
			StdLargeVec<Vecd> &pos_n_;
			ComplexShape &body_shape_;
			Kernel &kernel_;
			ParticleSpacingByBodyShape *particle_spacing_by_body_shape_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class RelaxationAccelerationComplex
		* @brief compute relaxation acceleration while consider the present of contact bodies
		* with considering contact interaction
		* this is usually used for fluid like bodies
		*/
		class RelaxationAccelerationComplex : public InteractionDynamics,
											  public RelaxDataDelegateComplex
		{
		public:
			explicit RelaxationAccelerationComplex(ComplexBodyRelation &body_complex_relation);
			virtual ~RelaxationAccelerationComplex(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &dvel_dt_;
			StdVec<StdLargeVec<Real> *> contact_Vol_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ShapeSurfaceBounding
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ShapeSurfaceBounding : public PartDynamicsByCell,
									 public RelaxDataDelegateSimple
		{
		public:
			ShapeSurfaceBounding(SPHBody &sph_body, NearShapeSurface &body_part);
			virtual ~ShapeSurfaceBounding(){};

		protected:
			StdLargeVec<Vecd> &pos_n_;
			LevelSetShape *level_set_shape_;
			Real constrained_distance_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ConstraintSurfaceParticles
		* @brief constrain surface particles by
		* map contrained particles to geometry face and
		* r = r + phi * norm (vector distance to face)
		*/
		class ConstraintSurfaceParticles : public PartSimpleDynamicsByParticle,
										   public RelaxDataDelegateSimple
		{
		public:
			ConstraintSurfaceParticles(SPHBody &sph_body, BodySurface &body_part);
			virtual ~ConstraintSurfaceParticles(){};

		protected:
			Real constrained_distance_;
			StdLargeVec<Vecd> &pos_n_;
			LevelSetShape *level_set_shape_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class RelaxationStepInner
		* @brief carry out particle relaxation step of particles within the body
		*/
		class RelaxationStepInner : public ParticleDynamics<void>
		{
		protected:
			RealBody *real_body_;
			BaseBodyRelationInner &inner_relation_;
			NearShapeSurface near_shape_surface_;

		public:
			explicit RelaxationStepInner(BaseBodyRelationInner &inner_relation, bool level_set_correction = false);
			virtual ~RelaxationStepInner(){};

			UniquePtr<RelaxationAccelerationInner> relaxation_acceleration_inner_;
			GetTimeStepSizeSquare get_time_step_square_;
			UpdateParticlePosition update_particle_position_;
			ShapeSurfaceBounding surface_bounding_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		* @class UpdateParticlePosition
		* @brief update the particle position for a time step
		*/
		class UpdateSolidParticlePosition : public ParticleDynamicsSimple, public solid_dynamics::SolidDataSimple
		{
		public:
			explicit UpdateSolidParticlePosition(SPHBody &sph_body);
			virtual ~UpdateSolidParticlePosition(){};

		protected:
			StdLargeVec<Vecd> &pos_0_, &pos_n_, &dvel_dt_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class SolidRelaxationStepInner
		* @brief carry out particle relaxation step of particles within the body
		*/
		class SolidRelaxationStepInner : public RelaxationStepInner, public solid_dynamics::SolidDataSimple
		{
		public:
			explicit SolidRelaxationStepInner(BaseBodyRelationInner &inner_relation, bool level_set_correction = false)
				: RelaxationStepInner(inner_relation, level_set_correction),
				  solid_dynamics::SolidDataSimple(*inner_relation.sph_body_),
				  update_solid_particle_position_(*real_body_){};
			virtual ~SolidRelaxationStepInner(){};

			UpdateSolidParticlePosition update_solid_particle_position_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};
	}
}
#endif //RELAX_DYNAMICS_H