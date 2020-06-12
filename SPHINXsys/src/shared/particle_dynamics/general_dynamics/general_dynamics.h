/**
* @file 	general_dynamics.h
* @brief 	This is the particle dynnamics apllicable for all type bodies
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "all_particle_dynamics.h"

namespace SPH
{
	/**
	* @class InitializeATimeStep
	* @brief initialize a time step for a body.
	* including initialize particle acceleration 
	* induced by viscous, gravity and other forces,
	* set number of ghost particles into zero.
	*/
	class InitializeATimeStep : public ParticleDynamicsSimple<SPHBody, BaseParticles>
	{
	protected:
		Gravity* gravity_;
		virtual void SetupDynamics(Real dt = 0.0) override;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
	public:
		InitializeATimeStep(SPHBody* body, Gravity* gravity = new Gravity(Vecd(0)));
		virtual ~InitializeATimeStep() {};
	};

	/**
	* @class RandomizePartilePosition
	* @brief Randomize the initial particle position
	*/
	class RandomizePartilePosition : public ParticleDynamicsSimple<SPHBody, BaseParticles>
	{
	protected:
		Real particle_spacing_;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
	public:
		RandomizePartilePosition(SPHBody* body);
		virtual ~RandomizePartilePosition() {};
	};

	/**
	* @class BoundingBodyDomain
	* @brief The base calss bounding particle position within a box body domain.
	*/
	class BoundingBodyDomain : public ParticleDynamicsByCells<SPHBody>
	{
		/** obtain the cells lower and upper boundy for the body domain. */
		void SetCellBounds();

	protected:
		/** lower and upper bound for checking. */
		Vecd body_lower_bound_, body_upper_bound_;
		Vecu body_lower_bound_cell_, body_upper_bound_cell_;

	public:
		BoundingBodyDomain(SPHBody* body);
		virtual ~BoundingBodyDomain() {};
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

		virtual void CheckLowerBound(size_t index_particle_i, Vecd& pnt, Real dt = 0.0);
		virtual void CheckUpperBound(size_t index_particle_i, Vecd& pnt, Real dt = 0.0);
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
		virtual void CheckLowerBound(size_t index_particle_i, Vecd& pnt, Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd& pnt, Real dt = 0.0) override;
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
			virtual void checkLowerBound(size_t index_particle_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_particle_i, Real dt = 0.0);
			InnerFunctor checking_bound_;
		public:
			Bounding(CellVector& bound_cells, SPHBody* body, int axis_direction, bool positive);
			virtual ~Bounding() {};
			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		class CreatingGhostParticles : public Bounding
		{
		protected:
			IndexVector& ghost_particles_;
			virtual void SetupDynamics(Real dt = 0.0) override { ghost_particles_.clear(); };
			virtual void checkLowerBound(size_t index_particle_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_particle_i, Real dt = 0.0) override;
		public:
			CreatingGhostParticles(IndexVector& ghost_particles, CellVector& bound_cells, 
				SPHBody* body, int axis_direction, bool positive);
			virtual ~CreatingGhostParticles() {};
			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};

		class UpdatingGhostStates : public BoundingInAxisDirection
		{
		protected:
			IndexVector& ghost_particles_;
			void updateForLowerBound(size_t index_particle_i, Real dt = 0.0);
			void updateForUpperBound(size_t index_particle_i, Real dt = 0.0);
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
	 * @brief  check whether paritcle velocity within a bound
	 */
	class VelocityBoundCheck 
		: public ParticleDynamicsReduce<bool, ReduceOR, SPHBody, BaseParticles>
	{
	protected:
		Real velocity_bound_;
		bool ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;

	public:
		VelocityBoundCheck(SPHBody* body, Real velocity_bound);
		virtual ~VelocityBoundCheck() {};
	};

	/**
	 * @class UpperFrontInXDirection
	 * @brief Get the upper front In X Direction for a SPH body
	 */
	class UpperFrontInXDirection : public ParticleDynamicsReduce<Real, ReduceMax, SPHBody>
	{
	protected:
		Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
	public:
		explicit UpperFrontInXDirection(SPHBody* body)
			: ParticleDynamicsReduce(body) {
			initial_reference_ = 0.0;
		};
		virtual ~UpperFrontInXDirection() {};
	};
}
