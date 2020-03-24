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
		Vecd initial_value_;
		virtual void SetupDynamics(Real dt = 0.0) override;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
	public:
		InitializeATimeStep(SPHBody* body);
		InitializeATimeStep(SPHBody* body, ExternalForce *external_force);
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
	* @class BoundingAxisDirection
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

		//cells in which particle checked for bounding
		CellVector lower_bound_cells_, upper_bound_cells_;

		virtual void CheckLowerBound(size_t index_particle_i, Vecd& pnt, Real dt = 0.0) = 0;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd& pnt, Real dt = 0.0) = 0;
	public:
		BoundingInAxisDirection(SPHBody* body, int axis_direction);
		virtual ~BoundingInAxisDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class PeriodicBoundingInAxisDirection
	* @brief Periodic bounding particle position in an axis direction
	*/
	class PeriodicBoundingInAxisDirection : public BoundingInAxisDirection
	{
	protected:
		Vecd periodic_translation_;

		virtual void CheckLowerBound(size_t index_particle_i, Vecd& pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd& pnt,
			Real dt = 0.0) override;
	public:
		PeriodicBoundingInAxisDirection(SPHBody* body, int axis_direction);
		virtual ~PeriodicBoundingInAxisDirection() {};
	};

	/**
	* @class PeriodicConditionInAxisDirection
	* @brief Periodic boundary condition in an axis direction
	*/
	class PeriodicConditionInAxisDirection
		: public PeriodicBoundingInAxisDirection
	{
	protected:
		virtual void CheckLowerBound(size_t index_particle_i, Vecd& pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd& pnt,
			Real dt = 0.0) override;
	public:

		PeriodicConditionInAxisDirection(SPHBody* body, int axis_direction)
			: PeriodicBoundingInAxisDirection(body, axis_direction) {};
		virtual ~PeriodicConditionInAxisDirection() {};

		/** This class is only implemented in sequential due to memory conflicts. */
		virtual void parallel_exec(Real dt = 0.0) override { exec(); };
	};

	/**
	* @class MirrorBoundingInAxisDirection
	* @brief Mirror bounding particle position and velocity in an axis direction
	*/
	class MirrorBoundingInAxisDirection : public BoundingBodyDomain
	{
	protected:
		/** the axis direction for bounding*/
		const int axis_;	
		/** the second axis according right hand rule. */
		const int second_axis_;
		/** the third axis according right hand rule. used only for 3d. */
		const int third_axis_;
		/** cells in which particle checked for bounding */
		CellVector bound_cells_;
		virtual void CheckLowerBound(size_t index_particle_i, Real dt = 0.0);
		virtual void CheckUpperBound(size_t index_particle_i, Real dt = 0.0);
		InnerFunctor checking_bound_;
	public:
		MirrorBoundingInAxisDirection(SPHBody* body, int axis_direction, bool positive);
		virtual ~MirrorBoundingInAxisDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class MirrorBoundaryConditionInAxisDirection
	* @brief Mirror boundary condition in an axis direction
	*/
	class MirrorBoundaryCreateGhostsInAxisDirection : public MirrorBoundingInAxisDirection
	{
	protected:
		IndexVector ghost_particles_;
		virtual void SetupDynamics(Real dt = 0.0) override { ghost_particles_.clear(); };
		virtual void CheckLowerBound(size_t index_particle_i, Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Real dt = 0.0) override;
	public:
		/**
		 * @brief Constructor.
		 * @param[in] fluid body.
		 * @param[in] body part by particles.
		 * @param[in] axis direction of in flow: 0, 1, 2 for x-, y- and z-axis.
		 * @param[in] direction sign of the inlfow: ture for positive direction.
		 */
		explicit MirrorBoundaryCreateGhostsInAxisDirection(FluidBody* body, int axis_direction, bool positive)
			: MirrorBoundingInAxisDirection(body, axis_direction, positive) {};
		virtual ~MirrorBoundaryCreateGhostsInAxisDirection() {};

		/** This class is only implemented in sequential due to memory conflicts. */
		virtual void parallel_exec(Real dt = 0.0) override { exec(); };
	};

	/**
	* @class MirrorBoundaryConditionInAxisDirection
	* @brief Mirror boundary condition in an axis direction
	*/
	class MirrorBoundaryConditionInAxisDirection 
		: public MirrorBoundaryCreateGhostsInAxisDirection
	{
	protected:
		void CheckLowerBoundUpdate(size_t index_particle_i, Real dt = 0.0);
		void CheckUpperBoundUpdate(size_t index_particle_i, Real dt = 0.0);
		InnerFunctor checking_bound_update_;
	public:
		/**
		 * @brief Constructor.
		 * @param[in] fluid body.
		 * @param[in] body part by particles.
		 * @param[in] axis direction of in flow: 0, 1, 2 for x-, y- and z-axis.
		 * @param[in] direction sign of the inlfow: ture for positive direction.
		 */
		explicit MirrorBoundaryConditionInAxisDirection(FluidBody* body, int axis_direction, bool positive);
		virtual ~MirrorBoundaryConditionInAxisDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
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
