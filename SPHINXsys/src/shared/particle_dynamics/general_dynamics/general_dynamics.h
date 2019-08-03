/**
* @file 	general_dynamics.h
* @brief 	This is the particle dynnamics apllicable for all type bodies
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "all_particle_dynamics.h"

namespace SPH
{
	/**
	* @class InitializeOtherAccelerations
	* @brief initialize particle acceleration
	*/
	class InitializeOtherAccelerations : public ParticleDynamicsSimple<SPHBody, Particles>
	{
	protected:
		Vecd initial_value_;
		virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
	public:
		InitializeOtherAccelerations(SPHBody* body);
		InitializeOtherAccelerations(SPHBody* body, ExternalForce *external_force);
		virtual ~InitializeOtherAccelerations() {};
	};

	/**
	* @class RandomizePartilePosition
	* @brief Randomize the initialize particle position
	*/
	class RandomizePartilePosition : public ParticleDynamicsSimple<SPHBody, Particles>
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
		//obtain the cells lower and upper boundy 
		//for the body domain
		void SetCellBounds();

	protected:
		//lower and upper bound for checking
		Vecd body_lower_bound_, body_upper_bound_;
		Vecu body_lower_bound_cell_, body_upper_bound_cell_;

		//dynamics of a particle
		//to be realized in specific algorithms
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;

	public:
		BoundingBodyDomain(SPHBody* body);
		virtual ~BoundingBodyDomain() {};
	};

	/**
	* @class BoundingInXDirection
	* @brief Bounding particle position in x direction, in genreal
	*/
	class BoundingInXDirection : public BoundingBodyDomain
	{
	protected:

		//cells in which particle checked for bounding
		StdVec<Vecu> lower_bound_cells_, upper_bound_cells_;

		//dynamics of a particle
		//to be realized in specific algorithms
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;

	public:
		BoundingInXDirection(SPHBody* body);
		virtual ~BoundingInXDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class BoundingInYDirection
	* @brief Bounding particle position in y direction, in genreal
	*/
	class BoundingInYDirection : public BoundingBodyDomain
	{
	protected:

		//cells in which particle checked for bounding
		StdVec<Vecu> lower_bound_cells_, upper_bound_cells_;

		//dynamics of a particle
		//to be realized in specific algorithms
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;

	public:
		BoundingInYDirection(SPHBody* body);
		virtual ~BoundingInYDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
	/**
	* @class BoundingInZDirection
	* @brief Bounding particle position in z direction, in genreal
	*/
	class BoundingInZDirection : public BoundingBodyDomain
	{
	protected:

		//cells in which particle checked for bounding
		StdVec<Vecu> lower_bound_cells_, upper_bound_cells_;

		//dynamics of a particle
		//to be realized in specific algorithms
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt = 0.0) = 0;

	public:
		BoundingInZDirection(SPHBody* body);
		virtual ~BoundingInZDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	* @class PeriodicBoundingInXDirection
	* @brief Periodic bounding particle position in x direction
	*/
	class PeriodicBoundingInXDirection : public BoundingInXDirection
	{
	protected:
		Vecd periodic_translation_;

		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
	public:

		PeriodicBoundingInXDirection(SPHBody* body);
		virtual ~PeriodicBoundingInXDirection() {};

	};

	/**
	* @class PeriodicBoundingInYDirection
	* @brief Periodic bounding particle position in y direction
	*/
	class PeriodicBoundingInYDirection : public BoundingInYDirection
	{
	protected:
		Vecd periodic_translation_;

		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
	public:

		PeriodicBoundingInYDirection(SPHBody* body);
		virtual ~PeriodicBoundingInYDirection() {};

	};

	/**
	* @class PeriodicBoundingInZDirection
	* @brief Periodic bounding particle position in z direction
	*/
	class PeriodicBoundingInZDirection : public BoundingInZDirection
	{
	protected:
		Vecd periodic_translation_;

		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
	public:

		PeriodicBoundingInZDirection(SPHBody* body);
		virtual ~PeriodicBoundingInZDirection() {};

	};

	/**
	* @class PeriodicConditionInXDirection
	* @brief Periodic boundary condition in x direction
	*/
	class PeriodicConditionInXDirection
		: public PeriodicBoundingInXDirection
	{
	protected:
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
	public:

		PeriodicConditionInXDirection(SPHBody* body);
		virtual ~PeriodicConditionInXDirection() {};

	};

	/**
	* @class PeriodicConditionInYDirection
	* @brief Periodic boundary condition in y direction
	*/
	class PeriodicConditionInYDirection
		: public PeriodicBoundingInYDirection
	{
	protected:
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
	public:

		PeriodicConditionInYDirection(SPHBody* body);
		virtual ~PeriodicConditionInYDirection() {};

	};

	/**
	* @class PeriodicConditionInZDirection
	* @brief Periodic boundary condition in z direction
	*/
	class PeriodicConditionInZDirection
		: public PeriodicBoundingInZDirection
	{
	protected:
		virtual void CheckLowerBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
		virtual void CheckUpperBound(size_t index_particle_i, Vecd pnt,
			Real dt = 0.0) override;
	public:

		PeriodicConditionInZDirection(SPHBody* body);
		virtual ~PeriodicConditionInZDirection() {};

	};

	/**
	 * @class VelocityBoundCheck
	 * @brief  check whether paritcle velocity within a bound
	 */
	class VelocityBoundCheck 
		: public ParticleDynamicsReduce<bool, ReduceOR, SPHBody, Particles>
	{
	protected:
		Real velocity_bound_;
		bool ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;

	public:
		VelocityBoundCheck(SPHBody* body, Real velocity_bound);
		virtual ~VelocityBoundCheck() {};
	};

	/**
	 * @class UpperBoundInXDirection
	 * @brief Get the Upper Bound In X Direction for a SPH body
	 */
	class UpperBoundInXDirection : public ParticleDynamicsReduce<Real, ReduceMax, SPHBody>
	{
	protected:
		Real ReduceFunction(size_t index_particle_i, Real dt = 0.0) override;
	public:
		explicit UpperBoundInXDirection(SPHBody* body)
			: ParticleDynamicsReduce(body) {
			initial_reference_ = 0.0;
		};
		virtual ~UpperBoundInXDirection() {};
	};
}
