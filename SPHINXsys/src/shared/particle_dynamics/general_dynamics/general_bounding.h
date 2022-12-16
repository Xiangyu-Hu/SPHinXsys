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
 * @file 	general_bounding.h
 * @brief 	This is the particle dynamics for domain bounding
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_BOUNDING_H
#define GENERAL_BOUNDING_H

#include "general_dynamics.h"

namespace SPH
{

	/**
	 * @class BoundingAlongAxis
	 * @brief Bounding particle position in along axis.
	 * The axis must be 0, 1 for 2d and 0, 1, 2 for 3d
	 */
	class BoundingAlongAxis : public BaseDynamics<void>,
							  public LocalDynamics,
							  public GeneralDataDelegateSimple
	{
	protected:
		const int axis_;			  /**< the axis directions for bounding*/
		BoundingBox bounding_bounds_; /**< lower and upper bound for checking. */
		StdLargeVec<Vecd> &pos_;
		BaseCellLinkedList &cell_linked_list_;
		Real cut_off_radius_max_; /**< maximum cut off radius to avoid boundary particle depletion */

	public:
		BoundingAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis);
		virtual ~BoundingAlongAxis(){};
	};

	/**
	 * @class BasePeriodicCondition
	 * @brief Base class for two different type periodic boundary conditions.
	 */
	class BasePeriodicCondition
	{
	protected:
		Vecd periodic_translation_;
		StdVec<CellLists> bound_cells_data_;
		Vecd setPeriodicTranslation(BoundingBox &bounding_bounds, int axis);

		/**
		 * @class PeriodicBounding
		 * @brief Periodic bounding particle position in an axis direction
		 */
		class PeriodicBounding : public BoundingAlongAxis
		{
		protected:
			Vecd &periodic_translation_;
			StdVec<CellLists> &bound_cells_data_;

			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);

		public:
			PeriodicBounding(Vecd &periodic_translation,
							 StdVec<CellLists> &bound_cells_data,
							 RealBody &real_body, BoundingBox bounding_bounds, int axis)
				: BoundingAlongAxis(real_body, bounding_bounds, axis),
				  periodic_translation_(periodic_translation),
				  bound_cells_data_(bound_cells_data){};
			virtual ~PeriodicBounding(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		BasePeriodicCondition(RealBody &real_body, BoundingBox bounding_bounds, int axis);
		virtual ~BasePeriodicCondition(){};
	};

	/**
	 * @class PeriodicConditionUsingCellLinkedList
	 * @brief The method imposing periodic boundary condition in an axis direction.
	 *	It includes two different steps, i.e. imposing periodic bounding and condition.
	 *	The first step is carried out before update cell linked list and
	 *	the second after the updating.
	 *	If the exec or parallel_exec is called directly, error message will be given.
	 */
	class PeriodicConditionUsingCellLinkedList : public BasePeriodicCondition
	{
	protected:
		/**
		 * @class PeriodicCellLinkedList
		 * @brief Periodic boundary condition in an axis direction
		 */
		class PeriodicCellLinkedList : public BoundingAlongAxis
		{
		protected:
			std::mutex mutex_cell_list_entry_; /**< mutex exclusion for memory conflict */
			Vecd &periodic_translation_;
			StdVec<CellLists> &bound_cells_data_;
			virtual void checkLowerBound(ListDataVector &cell_list_data, Real dt = 0.0);
			virtual void checkUpperBound(ListDataVector &cell_list_data, Real dt = 0.0);

		public:
			PeriodicCellLinkedList(Vecd &periodic_translation,
								   StdVec<CellLists> &bound_cells_data,
								   RealBody &real_body, BoundingBox bounding_bounds, int axis)
				: BoundingAlongAxis(real_body, bounding_bounds, axis),
				  periodic_translation_(periodic_translation),
				  bound_cells_data_(bound_cells_data){};
			;
			virtual ~PeriodicCellLinkedList(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		PeriodicConditionUsingCellLinkedList(RealBody &real_body, BoundingBox bounding_bounds, int axis)
			: BasePeriodicCondition(real_body, bounding_bounds, axis),
			  bounding_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis),
			  update_cell_linked_list_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis){};
		virtual ~PeriodicConditionUsingCellLinkedList(){};

		PeriodicBounding bounding_;
		PeriodicCellLinkedList update_cell_linked_list_;
	};

	/**
	 * @class PeriodicConditionUsingGhostParticles
	 * @brief The method imposing periodic boundary condition in an axis direction by using ghost particles.
	 *	It includes three different steps, i.e. imposing periodic bounding, creating ghosts and update ghost state.
	 *	The first step is carried out before update cell linked list and
	 *	the second and third after the updating.
	 *	If the exec or parallel_exec is called directly, error message will be given.
	 *  Note that, currently, this class is not for periodic condition in combined directions,
	 *  such as periodic condition in both x and y directions.
	 */
	class PeriodicConditionUsingGhostParticles : public BasePeriodicCondition
	{
	protected:
		StdVec<IndexVector> ghost_particles_;

		/**
		 * @class CreatPeriodicGhostParticles
		 * @brief create ghost particles in an axis direction
		 */
		class CreatPeriodicGhostParticles : public PeriodicBounding
		{
		protected:
			std::mutex mutex_create_ghost_particle_; /**< mutex exclusion for memory conflict */
			StdVec<IndexVector> &ghost_particles_;
			StdLargeVec<Real> &Vol_;
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;

		public:
			CreatPeriodicGhostParticles(Vecd &periodic_translation,
										StdVec<CellLists> &bound_cells_data,
										StdVec<IndexVector> &ghost_particles,
										RealBody &real_body, BoundingBox bounding_bounds, int axis)
				: PeriodicBounding(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis),
				  ghost_particles_(ghost_particles), Vol_(particles_->Vol_) {};
			virtual ~CreatPeriodicGhostParticles(){};
		};

		/**
		 * @class UpdatePeriodicGhostParticles
		 * @brief update ghost particles in an axis direction
		 */
		class UpdatePeriodicGhostParticles : public PeriodicBounding
		{
		protected:
			StdVec<IndexVector> &ghost_particles_;
			void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			void checkUpperBound(size_t index_i, Real dt = 0.0) override;

		public:
			UpdatePeriodicGhostParticles(Vecd &periodic_translation,
										 StdVec<CellLists> &bound_cells_data,
										 StdVec<IndexVector> &ghost_particles,
										 RealBody &real_body, BoundingBox bounding_bounds, int axis)
				: PeriodicBounding(periodic_translation, bound_cells_data, real_body, bounding_bounds, axis),
				  ghost_particles_(ghost_particles){};
			virtual ~UpdatePeriodicGhostParticles(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		PeriodicConditionUsingGhostParticles(RealBody &real_body, BoundingBox bounding_bounds, int axis)
			: BasePeriodicCondition(real_body, bounding_bounds, axis),
			  bounding_(periodic_translation_, bound_cells_data_, real_body, bounding_bounds, axis),
			  ghost_creation_(periodic_translation_, bound_cells_data_, ghost_particles_, real_body, bounding_bounds, axis),
			  ghost_update_(periodic_translation_, bound_cells_data_, ghost_particles_, real_body, bounding_bounds, axis)
		{
			ghost_particles_.resize(2);
		};

		virtual ~PeriodicConditionUsingGhostParticles(){};

		PeriodicBounding bounding_;
		CreatPeriodicGhostParticles ghost_creation_;
		UpdatePeriodicGhostParticles ghost_update_;
	};

	/**
	 * @class MirrorConditionAlongAxis
	 * @brief Mirror bounding particle position and velocity in an axis direction
	 *  Note that, currently, this class is not for mirror condition in combined directions,
	 *  such as mirror condition in both x and y directions.
	 */
	class MirrorConditionAlongAxis : public BoundingAlongAxis
	{
	protected:
		CellLists bound_cells_data_;
		IndexVector ghost_particles_;

		class MirrorBounding : public BoundingAlongAxis
		{
		protected:
			CellLists &bound_cells_data_; 
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);
			ParticleFunctor checking_bound_;

			StdLargeVec<Vecd> &vel_;
			void mirrorAlongAxis(size_t particle_index_i, Vecd body_bound, int axis);

		public:
			MirrorBounding(CellLists &bound_cells_data, RealBody &real_body,
						   BoundingBox bounding_bounds, int axis, bool positive);
			virtual ~MirrorBounding(){};
			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		 * @class CreatingMirrorGhostParticles
		 * @brief ghost particle created according to its corresponding real particle
		 */
		class CreatingMirrorGhostParticles : public MirrorBounding
		{
		protected:
			IndexVector &ghost_particles_;
			virtual void setupDynamics(Real dt = 0.0) override { ghost_particles_.clear(); };
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;

		public:
			CreatingMirrorGhostParticles(IndexVector &ghost_particles,
										 CellLists &bound_cells_data, RealBody &real_body,
										 BoundingBox bounding_bounds, int axis, bool positive);
			virtual ~CreatingMirrorGhostParticles(){};
			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};

		/**
		 * @class UpdatingMirrorGhostStates
		 * @brief the state of a ghost particle updated according to its corresponding real particle
		 */
		class UpdatingMirrorGhostStates : public MirrorBounding
		{
		protected:
			IndexVector &ghost_particles_;
			void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			void checkUpperBound(size_t index_i, Real dt = 0.0) override;
			ParticleFunctor checking_bound_update_;

		public:
			UpdatingMirrorGhostStates(IndexVector &ghost_particles,
									  CellLists &bound_cells_data,
									  RealBody &real_body, BoundingBox bounding_bounds, int axis, bool positive);
			virtual ~UpdatingMirrorGhostStates(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		MirrorBounding bounding_;
		CreatingMirrorGhostParticles creating_ghost_particles_;
		UpdatingMirrorGhostStates updating_ghost_states_;

		MirrorConditionAlongAxis(RealBody &real_body, BoundingBox bounding_bounds, int axis, bool positive);
		virtual ~MirrorConditionAlongAxis(){};

		virtual void exec(Real dt = 0.0) override{};
		virtual void parallel_exec(Real dt = 0.0) override{};
	};
}
#endif // GENERAL_BOUNDING_H
