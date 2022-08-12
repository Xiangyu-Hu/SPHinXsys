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

#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "cell_linked_list.h"
#include "all_body_relations.h"
#include "general_dynamics.h"

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
			StdLargeVec<Vecd> &acc_;
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
			StdLargeVec<Vecd> &acc_;
			StdLargeVec<Vecd> &pos_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class RelaxationAccelerationInnerWithLevelSetCorrection
		 * @brief we constrain particles to a level function representing the interafce.
		 */
		class RelaxationAccelerationInnerWithLevelSetCorrection : public RelaxationAccelerationInner
		{
		public:
			explicit RelaxationAccelerationInnerWithLevelSetCorrection(
				BaseBodyRelationInner &inner_relation);
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
			StdLargeVec<Vecd> &pos_, &acc_;
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
			StdLargeVec<Vecd> &pos_;
			Shape &body_shape_;
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
			StdLargeVec<Vecd> &acc_, &pos_;
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
			StdLargeVec<Vecd> &pos_;
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
			StdLargeVec<Vecd> &pos_;
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
			explicit RelaxationStepInner(BaseBodyRelationInner &inner_relation,
										 bool level_set_correction = false);
			virtual ~RelaxationStepInner(){};

			UniquePtr<RelaxationAccelerationInner> relaxation_acceleration_inner_;
			GetTimeStepSizeSquare get_time_step_square_;
			UpdateParticlePosition update_particle_position_;
			ShapeSurfaceBounding surface_bounding_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		 * @class RelaxationAccelerationComplexWithLevelSetCorrection
		 * @brief compute relaxation acceleration while consider the present of contact bodies
		 * with considering contact interaction
		 * this is usually used for fluid like bodies
		 * we constrain particles with a level-set correction function when the fluid boundary is not contacted with solid.
		 */
		class RelaxationAccelerationComplexWithLevelSetCorrection : public RelaxationAccelerationComplex
		{
		public:
			RelaxationAccelerationComplexWithLevelSetCorrection(
				ComplexBodyRelation &body_complex_relation, const std::string &shape_name);
			virtual ~RelaxationAccelerationComplexWithLevelSetCorrection(){};

		protected:
			LevelSetShape *level_set_shape_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class RelaxationStepComplex
		 * @brief carry out particle relaxation step of particles within multi bodies
		 */
		class RelaxationStepComplex : public ParticleDynamics<void>
		{
		protected:
			RealBody *real_body_;
			ComplexBodyRelation &complex_relation_;
			NearShapeSurface near_shape_surface_;

		public:
			explicit RelaxationStepComplex(ComplexBodyRelation &body_complex_relation,
										   const std::string &shape_name, bool level_set_correction = false);
			virtual ~RelaxationStepComplex(){};

			UniquePtr<RelaxationAccelerationComplex> relaxation_acceleration_complex_;
			GetTimeStepSizeSquare get_time_step_square_;
			UpdateParticlePosition update_particle_position_;
			ShapeSurfaceBounding surface_bounding_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		* @class ShellMidSurfaceBounding
		* @brief constrain particles by constraining particles to mid-surface.
		* Note that level_set_refinement_ratio should be smaller than particle_spacing_ref_ / (0.05 * thickness_)
		* because if level_set_refinement_ratio > particle_spacing_ref_ / (0.05 * thickness_), 
		* there will be no level set field.
		*/
		class ShellMidSurfaceBounding : public PartDynamicsByCell,
										public RelaxDataDelegateInner
		{
		public:
			ShellMidSurfaceBounding(SPHBody &body, NearShapeSurface &body_part, BaseBodyRelationInner &inner_relation,
									Real thickness, Real level_set_refinement_ratio);
			virtual ~ShellMidSurfaceBounding(){};

		protected:
			StdLargeVec<Vecd> &pos_;
			Real constrained_distance_;
			LevelSetShape *level_set_shape_;
			Real particle_spacing_ref_, thickness_, level_set_refinement_ratio_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ShellNormalDirectionPrediction
		 * @brief prodict the normal direction of shell particles.
		 */
		class ShellNormalDirectionPrediction : public ParticleDynamics<void>
		{
			const Real convergence_criterion_;
			const Real consistency_criterion_;

			void predictNormalDirection();
			void correctNormalDirection();

		public:
			explicit ShellNormalDirectionPrediction(BaseBodyRelationInner &inner_relation, 
				Real thickness, Real consistency_criterion = cos(Pi / 20.0));
			virtual ~ShellNormalDirectionPrediction(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };

		protected:
			class NormalPrediction : public RelaxDataDelegateSimple
			{
				Real thickness_;
				LevelSetShape *level_set_shape_;
				StdLargeVec<Vecd> &pos_, &n_, n_temp_;

			public:
				NormalPrediction(SPHBody &sph_body, Real thickness);
				virtual ~NormalPrediction(){};
				void update(size_t index_i, Real dt = 0.0);
			};

			class PredictionConvergenceCheck : public ParticleDynamicsReduce<bool, ReduceAND>,
											   public RelaxDataDelegateSimple
			{
			public:
				PredictionConvergenceCheck(SPHBody &sph_body, Real convergence_criterion);
				virtual ~PredictionConvergenceCheck(){};

			protected:
				const Real convergence_criterion_;
				StdLargeVec<Vecd> &n_, &n_temp_;
				bool ReduceFunction(size_t index_i, Real dt = 0.0) override;
			};

			class ConsistencyCorrection : public InteractionDynamics, public RelaxDataDelegateInner
			{
			public:
				explicit ConsistencyCorrection(BaseBodyRelationInner &inner_relation, Real consistency_criterion);
				virtual ~ConsistencyCorrection(){};

				/** only implement sequential version now. */
				virtual void parallel_exec(Real dt = 0.0) override { exec(); };

			protected:
				const Real consistency_criterion_;
				StdLargeVec<int> updated_indicator_; /**> 0 not updated, 1 updated with reliable prediction, 2 updated from a reliable neighbor */
				StdLargeVec<Vecd> &n_;
				virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			};

			class ConsistencyUpdatedCheck : public ParticleDynamicsReduce<bool, ReduceAND>,
								 public RelaxDataDelegateSimple
			{
			public:
				explicit ConsistencyUpdatedCheck(SPHBody &sph_body);
				virtual ~ConsistencyUpdatedCheck(){};

			protected:
				StdLargeVec<int> &updated_indicator_;
				bool ReduceFunction(size_t index_i, Real dt = 0.0) override;
			};

			class SmoothingNormal : public ParticleSmoothing<Vecd>
			{
			public:
				explicit SmoothingNormal(BaseBodyRelationInner &inner_relation);
				virtual ~SmoothingNormal(){};

			protected:
				virtual void Update(size_t index_i, Real dt = 0.0) override;
			};

			SimpleDynamics<NormalPrediction> normal_prediction_;
			PredictionConvergenceCheck normal_prediction_convergence_check_;
			ConsistencyCorrection consistency_correction_;
			ConsistencyUpdatedCheck consistency_updated_check_;
			SmoothingNormal smoothing_normal_;
		};

		/**
		 * @class ShellRelaxationStepInner
		 * @brief carry out particle relaxation step of particles within the shell body
		 */
		class ShellRelaxationStepInner : public RelaxationStepInner
		{
		public:
			explicit ShellRelaxationStepInner(BaseBodyRelationInner &inner_relation, Real thickness,
											  Real level_set_refinement_ratio, bool level_set_correction = false);
			virtual ~ShellRelaxationStepInner(){};

			UpdateParticlePosition update_shell_particle_position_;
			ShellMidSurfaceBounding mid_surface_bounding_;

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};
	}
}
#endif // RELAX_DYNAMICS_H