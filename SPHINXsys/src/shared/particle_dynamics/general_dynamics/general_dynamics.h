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
* @file 	general_dynamics.h
* @brief 	This is the particle dynamics aplliable for all type bodies
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "all_particle_dynamics.h"

#include <limits>

namespace SPH
{
	typedef DataDelegateSimple<SPHBody, BaseParticles> GeneralDataDelegateSimple;

	/**
	* @class InitializeATimeStep
	* @brief initialize a time step for a body.
	* including initialize particle acceleration 
	* induced by viscous, gravity and other forces,
	* set number of ghost particles into zero.
	*/
	class InitializeATimeStep 
		: public ParticleDynamicsSimple, public GeneralDataDelegateSimple
	{
	public:
		InitializeATimeStep(SPHBody* body, Gravity* gravity = new Gravity(Vecd(0)));
		virtual ~InitializeATimeStep() {};
	protected:
		StdLargeVec<Vecd>& pos_n_,& dvel_dt_others_;
		Gravity* gravity_;
		virtual void setupDynamics(Real dt = 0.0) override;
		virtual void Update(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class RandomizePartilePosition
	* @brief Randomize the initial particle position
	*/
	class RandomizePartilePosition
		: public ParticleDynamicsSimple, public GeneralDataDelegateSimple
	{
	public:
		RandomizePartilePosition(SPHBody* body);
		virtual ~RandomizePartilePosition() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		Real particle_spacing_;
		virtual void Update(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class BoundingBodyDomain
	* @brief The base class bounding particle position within a box body domain.
	*/
	class BoundingBodyDomain
		: public ParticleDynamics<void>, public GeneralDataDelegateSimple
	{
	public:
		BoundingBodyDomain(SPHBody* body);
		virtual ~BoundingBodyDomain() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		matrix_cell cell_linked_lists_;
		Vecu number_of_cells_;
		Real cell_spacing_;
		Vecd mesh_lower_bound_;

		/** lower and upper bound for checking. */
		Vecd body_lower_bound_, body_upper_bound_;
		Vecu body_lower_bound_cell_, body_upper_bound_cell_;
	private:
		/** obtain the cells lower and upper boundary for the body domain. */
		void SetCellBounds();
	};

	/**
	* @class BoundingInAxisDirection
	* @brief Bounding particle position in a axis direction.
	* The axis_direction must be 0, 1 for 2d and 0, 1, 2 for 3d
	*/
	class BoundingInAxisDirection : public BoundingBodyDomain
	{
	protected:
		/** the axis direction for bounding*/
		const int axis_;
		/** the second axis according right hand rule. */
		const int second_axis_;
		/** the third axis according right hand rule. used only for 3d. */
		const int third_axis_;
	public:
		BoundingInAxisDirection(SPHBody* body, int axis_direction);
		virtual ~BoundingInAxisDirection() {};
	};

	/**
	* @class PeriodicBoundingInAxisDirection
	* @brief Periodic bounding particle position in an axis direction
	*/
	class PeriodicBoundingInAxisDirection : public BoundingInAxisDirection
	{
	protected:
		Vecd periodic_translation_;
		//cells in which particle checked for bounding
		CellVector lower_bound_cells_, upper_bound_cells_;

		/**compute the distance for periodic translation. */
		void setPeriodicTranslation();

		virtual void CheckLowerBound(size_t index_i, Vecd& pnt, Real dt = 0.0);
		virtual void CheckUpperBound(size_t index_i, Vecd& pnt, Real dt = 0.0);
	public:
		PeriodicBoundingInAxisDirection(SPHBody* body, int axis_direction);
		virtual ~PeriodicBoundingInAxisDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class PeriodicConditionInAxisDirection
	* @brief Periodic boundary condition in an axis direction
	*/
	class PeriodicConditionInAxisDirection
		: public PeriodicBoundingInAxisDirection
	{
	protected:
		virtual void CheckLowerBound(size_t index_i, Vecd& pnt, Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_i, Vecd& pnt, Real dt = 0.0) override;
	public:

		PeriodicConditionInAxisDirection(SPHBody* body, int axis_direction)
			: PeriodicBoundingInAxisDirection(body, axis_direction) {};
		virtual ~PeriodicConditionInAxisDirection() {};

		/** This class is only implemented in sequential due to memory conflicts. */
		virtual void parallel_exec(Real dt = 0.0) override { exec(); };
	};

	/**
	* @class MirrorBoundaryConditionInAxisDirection
	* @brief Mirror bounding particle position and velocity in an axis direction
	*/
	class MirrorBoundaryConditionInAxisDirection : public BoundingInAxisDirection
	{
	protected:
		/** cells in which particle checked for bounding */
		CellVector bound_cells_;
		/** ghost particles createded for impose boundary condition. */
		IndexVector ghost_particles_;

		class Bounding : public BoundingInAxisDirection
		{
		protected:
			CellVector& bound_cells_;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);
			InnerFunctor checking_bound_;
		public:
			Bounding(CellVector& bound_cells, SPHBody* body, int axis_direction, bool positive);
			virtual ~Bounding() {};
			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		* @class CreatingGhostParticles
		* @brief ghost particle created according to its corresponding real particle
		*/
		class CreatingGhostParticles : public Bounding
		{
		protected:
			IndexVector& ghost_particles_;
			virtual void setupDynamics(Real dt = 0.0) override { ghost_particles_.clear(); };
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;
		public:
			CreatingGhostParticles(IndexVector& ghost_particles, CellVector& bound_cells, 
				SPHBody* body, int axis_direction, bool positive);
			virtual ~CreatingGhostParticles() {};
			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};

		/**
		* @class UpdatingGhostStates
		* @brief the state of a ghost particle updated according to its corresponding real particle
		*/
		class UpdatingGhostStates : public BoundingInAxisDirection
		{
		protected:
			IndexVector& ghost_particles_;
			void updateForLowerBound(size_t index_i, Real dt = 0.0);
			void updateForUpperBound(size_t index_i, Real dt = 0.0);
			InnerFunctor checking_bound_update_;
		public:
			UpdatingGhostStates(IndexVector& ghost_particles,
				SPHBody* body, int axis_direction, bool positive);
			virtual ~UpdatingGhostStates() {};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		Bounding bounding_;
		CreatingGhostParticles creating_ghost_particles_;
		UpdatingGhostStates updating_ghost_states_;

		MirrorBoundaryConditionInAxisDirection(SPHBody* body, int axis_direction, bool positive);
		virtual ~MirrorBoundaryConditionInAxisDirection() {};

		virtual void exec(Real dt = 0.0) override {};
		virtual void parallel_exec(Real dt = 0.0) override {};
	};

	/**
	 * @class VelocityBoundCheck
	 * @brief  check whether particle velocity within a given bound
	 */
	class VelocityBoundCheck : 
		public ParticleDynamicsReduce<bool, ReduceOR>, 
		public GeneralDataDelegateSimple
	{
	public:
		VelocityBoundCheck(SPHBody* body, Real velocity_bound);
		virtual ~VelocityBoundCheck() {};
	protected:
		StdLargeVec<Vecd>& vel_n_;
		Real velocity_bound_;
		bool ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class UpperFrontInXDirection
	 * @brief Get the upper front In X Direction for a SPH body
	 */
	class UpperFrontInXDirection : 
		public ParticleDynamicsReduce<Real, ReduceMax>, 
		public GeneralDataDelegateSimple
	{
	public:
		explicit UpperFrontInXDirection(SPHBody* body);
		virtual ~UpperFrontInXDirection() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class MaximumSpeed
	 * @brief Get the maximum particle speed in a SPH body
	 */
	class MaximumSpeed : 
		public ParticleDynamicsReduce<Real, ReduceMax>, 
		public GeneralDataDelegateSimple
	{
	public:
		explicit MaximumSpeed(SPHBody* body);
		virtual ~MaximumSpeed() {};
	protected:
		StdLargeVec<Vecd>& vel_n_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class BodyLowerBound
	* @brief the lower bound of a body by reduced particle positions.
	*/
	class BodyLowerBound : 
		public  ParticleDynamicsReduce<Vecd, ReduceLowerBound>,
		public GeneralDataDelegateSimple
	{
	public:
		explicit BodyLowerBound(SPHBody* body);
		virtual ~BodyLowerBound() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class BodyUpperBound
	 * @brief the upper bound of a body by reduced particle positions.
	 */
	class BodyUpperBound : 
		public  ParticleDynamicsReduce<Vecd, ReduceUpperBound>,
		public GeneralDataDelegateSimple
	{
	public:
		explicit BodyUpperBound(SPHBody* body);
		virtual ~BodyUpperBound() {};
	protected:
		StdLargeVec<Vecd>& pos_n_;
		Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};
}
