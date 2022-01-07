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
*/

#ifndef GENERAL_DYNAMICS_H
#define GENERAL_DYNAMICS_H

#include "all_particle_dynamics.h"

#include <limits>

namespace SPH
{
	typedef DataDelegateSimple<SPHBody, BaseParticles> GeneralDataDelegateSimple;
	typedef DataDelegateContact<SPHBody, BaseParticles, BaseMaterial,
								SPHBody, BaseParticles, BaseMaterial, DataDelegateEmptyBase>
		GeneralDataDelegateContact;
	/**
	* @class TimeStepInitialization
	* @brief initialize a time step for a body.
	* including initialize particle acceleration 
	* induced by viscous, gravity and other forces,
	* set the number of ghost particles into zero.
	*/
	class TimeStepInitialization
		: public ParticleDynamicsSimple,
		  public GeneralDataDelegateSimple
	{
	private:
		UniquePtrKeeper<Gravity> gravity_ptr_keeper_;

	public:
		explicit TimeStepInitialization(SPHBody &sph_body);
		TimeStepInitialization(SPHBody &sph_body, Gravity &gravity);
		virtual ~TimeStepInitialization(){};

	protected:
		StdLargeVec<Vecd> &pos_n_, &dvel_dt_prior_;
		Gravity *gravity_;
		virtual void setupDynamics(Real dt = 0.0) override;
		virtual void Update(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class RandomizePartilePosition
	* @brief Randomize the initial particle position
	*/
	class RandomizePartilePosition
		: public ParticleDynamicsSimple,
		  public GeneralDataDelegateSimple
	{
	public:
		RandomizePartilePosition(SPHBody &sph_body);
		virtual ~RandomizePartilePosition(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Real randomize_scale_;
		virtual void Update(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class BoundingInAxisDirection
	* @brief Bounding particle position in a axis direction.
	* The axis_direction must be 0, 1 for 2d and 0, 1, 2 for 3d
	*/
	class BoundingInAxisDirection : public ParticleDynamics<void>, public GeneralDataDelegateSimple
	{
	protected:
		const int axis_;				 /**< the axis directions for bounding*/
		BoundingBox body_domain_bounds_; /**< lower and upper bound for checking. */
		StdLargeVec<Vecd> &pos_n_;
		BaseCellLinkedList *cell_linked_list_;
		Real cut_off_radius_max_; /**< maximum cut off radius to avoid boundary particle depletion */
	public:
		BoundingInAxisDirection(RealBody &real_body, int axis_direction);
		virtual ~BoundingInAxisDirection(){};
	};

	/**
	 * @class PeriodicConditionInAxisDirection
	 * @brief Base class for two different type periodic boundary conditions.
	 */
	class PeriodicConditionInAxisDirection
	{
	protected:
		Vecd periodic_translation_;
		StdVec<CellLists> bound_cells_;
		void setPeriodicTranslation(BoundingBox &body_domain_bounds, int axis_direction);

		/**
		* @class PeriodicBounding
		* @brief Periodic bounding particle position in an axis direction
		*/
		class PeriodicBounding : public BoundingInAxisDirection
		{
		protected:
			Vecd &periodic_translation_;
			StdVec<CellLists> &bound_cells_;

			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);

		public:
			PeriodicBounding(Vecd &periodic_translation,
							 StdVec<CellLists> &bound_cells, RealBody &real_body, int axis_direction)
				: BoundingInAxisDirection(real_body, axis_direction), periodic_translation_(periodic_translation),
				  bound_cells_(bound_cells){};
			virtual ~PeriodicBounding(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		* @class PeriodicCondition
		* @brief implement periodic condition in an axis direction
		*/
		class PeriodicCondition : public BoundingInAxisDirection
		{
		protected:
			Vecd &periodic_translation_;
			StdVec<CellLists> &bound_cells_;

			virtual void checkLowerBound(ListData &list_data, Real dt = 0.0) = 0;
			virtual void checkUpperBound(ListData &list_data, Real dt = 0.0) = 0;

		public:
			PeriodicCondition(Vecd &periodic_translation,
							  StdVec<CellLists> &bound_cells, RealBody &real_body, int axis_direction)
				: BoundingInAxisDirection(real_body, axis_direction), periodic_translation_(periodic_translation),
				  bound_cells_(bound_cells){};
			virtual ~PeriodicCondition(){};

			/** This class is only implemented in sequential due to memory conflicts.
			 * Because the cell list data is not concurrent vector.
			 */
			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};

	public:
		PeriodicConditionInAxisDirection(RealBody &real_body, int axis_direction);
		virtual ~PeriodicConditionInAxisDirection(){};
	};

	/**
	 * @class PeriodicConditionInAxisDirectionUsingCellLinkedList
	 * @brief The method imposing periodic boundary condition in an axis direction.
	 *	It includes two different steps, i.e. imposing periodic bounding and condition.
	 *	The first step is carried out before update cell linked list and
	 *	the second after the updating.
	 *	If the exec or parallel_exec is called directly, error message will be given.
	 */
	class PeriodicConditionInAxisDirectionUsingCellLinkedList : public PeriodicConditionInAxisDirection
	{
	protected:
		/**
		* @class PeriodicCondition
		* @brief Periodic boundary condition in an axis direction
		*/
		class PeriodicCellLinkedList : public PeriodicCondition
		{
		protected:
			virtual void checkLowerBound(ListData &list_data, Real dt = 0.0) override;
			virtual void checkUpperBound(ListData &list_data, Real dt = 0.0) override;

		public:
			PeriodicCellLinkedList(Vecd &periodic_translation,
								   StdVec<CellLists> &bound_cells, RealBody &real_body, int axis_direction)
				: PeriodicCondition(periodic_translation, bound_cells, real_body, axis_direction){};
			virtual ~PeriodicCellLinkedList(){};
		};

	public:
		PeriodicConditionInAxisDirectionUsingCellLinkedList(RealBody &real_body, int axis_direction)
			: PeriodicConditionInAxisDirection(real_body, axis_direction),
			  bounding_(this->periodic_translation_, this->bound_cells_, real_body, axis_direction),
			  update_cell_linked_list_(this->periodic_translation_, this->bound_cells_, real_body, axis_direction){};
		virtual ~PeriodicConditionInAxisDirectionUsingCellLinkedList(){};

		PeriodicBounding bounding_;
		PeriodicCellLinkedList update_cell_linked_list_;
	};

	/**
    * @class OpenBoundaryConditionInAxisDirection
	* @brief In open boundary case, we transfer fluid particles to buffer particles at outlet
	* @brief int axis_direction is used to choose direction in coordinate
	* @brief bool positive is used to choose upper or lower bound in your chosen direction
	*/
	class OpenBoundaryConditionInAxisDirection
	{
	protected:
		StdVec<CellLists> bound_cells_;

		class ParticleTypeTransfer : public BoundingInAxisDirection
		{
		protected:
			StdVec<CellLists> &bound_cells_;
			ParticleFunctor checking_bound_;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);

		public:
			ParticleTypeTransfer(StdVec<CellLists> &bound_cells, RealBody &real_body, int axis_direction, bool positive)
				: BoundingInAxisDirection(real_body, axis_direction),
				  bound_cells_(bound_cells)
			{
				checking_bound_ = positive ? std::bind(&OpenBoundaryConditionInAxisDirection::ParticleTypeTransfer::checkUpperBound, this, _1, _2)
										   : std::bind(&OpenBoundaryConditionInAxisDirection::ParticleTypeTransfer::checkLowerBound, this, _1, _2);
			};
			virtual ~ParticleTypeTransfer(){};

			/** This class is only implemented in sequential due to memory conflicts.
			* Because the cell list data is not concurrent vector.
			*/
			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};

	public:
		OpenBoundaryConditionInAxisDirection(RealBody &real_body, int axis_direction, bool positive);
		virtual ~OpenBoundaryConditionInAxisDirection(){};

		ParticleTypeTransfer particle_type_transfer;
	};

	/**
	 * @class PeriodicConditionInAxisDirectionUsingGhostParticles
	 * @brief The method imposing periodic boundary condition in an axis direction by using ghost particles.
	 *	It includes three different steps, i.e. imposing periodic bounding, creating ghosts and update ghost state.
	 *	The first step is carried out before update cell linked list and
	 *	the second and third after the updating.
	 *	If the exec or parallel_exec is called directly, error message will be given.
	 *  Note that, currently, this class is not for periodic condition in combined directions,
	 *  such as periodic condition in both x and y directions.
	 */
	class PeriodicConditionInAxisDirectionUsingGhostParticles : public PeriodicConditionInAxisDirection
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
			StdVec<IndexVector> &ghost_particles_;
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;

		public:
			CreatPeriodicGhostParticles(Vecd &periodic_translation, StdVec<CellLists> &bound_cells,
										StdVec<IndexVector> &ghost_particles, RealBody &real_body, int axis_direction)
				: PeriodicBounding(periodic_translation, bound_cells, real_body, axis_direction),
				  ghost_particles_(ghost_particles){};
			virtual ~CreatPeriodicGhostParticles(){};

			/** This class is only implemented in sequential due to memory conflicts.
			 * Because creating ghost particle allocate memory.
			 */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
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
			UpdatePeriodicGhostParticles(Vecd &periodic_translation, StdVec<CellLists> &bound_cells,
										 StdVec<IndexVector> &ghost_particles, RealBody &real_body, int axis_direction)
				: PeriodicBounding(periodic_translation, bound_cells, real_body, axis_direction),
				  ghost_particles_(ghost_particles){};
			virtual ~UpdatePeriodicGhostParticles(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		PeriodicConditionInAxisDirectionUsingGhostParticles(RealBody &real_body, int axis_direction)
			: PeriodicConditionInAxisDirection(real_body, axis_direction),
			  bounding_(this->periodic_translation_, this->bound_cells_, real_body, axis_direction),
			  ghost_creation_(this->periodic_translation_, this->bound_cells_, this->ghost_particles_, real_body, axis_direction),
			  ghost_update_(this->periodic_translation_, this->bound_cells_, this->ghost_particles_, real_body, axis_direction)
		{
			ghost_particles_.resize(2);
		};

		virtual ~PeriodicConditionInAxisDirectionUsingGhostParticles(){};

		PeriodicBounding bounding_;
		CreatPeriodicGhostParticles ghost_creation_;
		UpdatePeriodicGhostParticles ghost_update_;
	};

	/**
	* @class MirrorBoundaryConditionInAxisDirection
	* @brief Mirror bounding particle position and velocity in an axis direction
	*  Note that, currently, this class is not for mirror condition in combined directions,
	*  such as mirror condition in both x and y directions.
	*/
	class MirrorBoundaryConditionInAxisDirection : public BoundingInAxisDirection
	{
	protected:
		CellLists bound_cells_;
		IndexVector ghost_particles_;

		class MirrorBounding : public BoundingInAxisDirection
		{
		protected:
			CellLists &bound_cells_;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);
			ParticleFunctor checking_bound_;

			StdLargeVec<Vecd> &vel_n_;
			void mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction);

		public:
			MirrorBounding(CellLists &bound_cells, RealBody &real_body, int axis_direction, bool positive);
			virtual ~MirrorBounding(){};
			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		* @class CreatingGhostParticles
		* @brief ghost particle created according to its corresponding real particle
		*/
		class CreatingGhostParticles : public MirrorBounding
		{
		protected:
			IndexVector &ghost_particles_;
			virtual void setupDynamics(Real dt = 0.0) override { ghost_particles_.clear(); };
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;

		public:
			CreatingGhostParticles(IndexVector &ghost_particles, CellLists &bound_cells,
								   RealBody &real_body, int axis_direction, bool positive);
			virtual ~CreatingGhostParticles(){};
			/** This class is only implemented in sequential due to memory conflicts. */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};

		/**
		* @class UpdatingGhostStates
		* @brief the state of a ghost particle updated according to its corresponding real particle
		*/
		class UpdatingGhostStates : public MirrorBounding
		{
		protected:
			IndexVector &ghost_particles_;
			void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			void checkUpperBound(size_t index_i, Real dt = 0.0) override;
			ParticleFunctor checking_bound_update_;

		public:
			UpdatingGhostStates(IndexVector &ghost_particles, CellLists &bound_cells,
								RealBody &real_body, int axis_direction, bool positive);
			virtual ~UpdatingGhostStates(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		MirrorBounding bounding_;
		CreatingGhostParticles creating_ghost_particles_;
		UpdatingGhostStates updating_ghost_states_;

		MirrorBoundaryConditionInAxisDirection(RealBody &real_body, int axis_direction, bool positive);
		virtual ~MirrorBoundaryConditionInAxisDirection(){};

		virtual void exec(Real dt = 0.0) override{};
		virtual void parallel_exec(Real dt = 0.0) override{};
	};

	/**
	 * @class VelocityBoundCheck
	 * @brief  check whether particle velocity within a given bound
	 */
	class VelocityBoundCheck : public ParticleDynamicsReduce<bool, ReduceOR>,
							   public GeneralDataDelegateSimple
	{
	public:
		VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound);
		virtual ~VelocityBoundCheck(){};

	protected:
		StdLargeVec<Vecd> &vel_n_;
		Real velocity_bound_;
		bool ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class UpperFrontInXDirection
	 * @brief Get the upper front In X Direction for a SPH body
	 */
	class UpperFrontInXDirection : public ParticleDynamicsReduce<Real, ReduceMax>,
								   public GeneralDataDelegateSimple
	{
	public:
		explicit UpperFrontInXDirection(SPHBody &sph_body);
		virtual ~UpperFrontInXDirection(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class MaximumSpeed
	 * @brief Get the maximum particle speed in a SPH body
	 */
	class MaximumSpeed : public ParticleDynamicsReduce<Real, ReduceMax>,
						 public GeneralDataDelegateSimple
	{
	public:
		explicit MaximumSpeed(SPHBody &sph_body);
		virtual ~MaximumSpeed(){};

	protected:
		StdLargeVec<Vecd> &vel_n_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class BodyLowerBound
	* @brief the lower bound of a body by reduced particle positions.
	*/
	class BodyLowerBound : public ParticleDynamicsReduce<Vecd, ReduceLowerBound>,
						   public GeneralDataDelegateSimple
	{
	public:
		explicit BodyLowerBound(SPHBody &sph_body);
		virtual ~BodyLowerBound(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	 * @class BodyUpperBound
	 * @brief the upper bound of a body by reduced particle positions.
	 */
	class BodyUpperBound : public ParticleDynamicsReduce<Vecd, ReduceUpperBound>,
						   public GeneralDataDelegateSimple
	{
	public:
		explicit BodyUpperBound(SPHBody &sph_body);
		virtual ~BodyUpperBound(){};

	protected:
		StdLargeVec<Vecd> &pos_n_;
		Vecd ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};

	/**
	* @class BodySummation
	* @brief Compute the summation of  a particle variable in a body
	*/
	template <int DataTypeIndex, typename VariableType>
	class BodySummation : public ParticleDynamicsReduce<VariableType, ReduceSum<VariableType>>,
						  public GeneralDataDelegateSimple
	{
	public:
		explicit BodySummation(SPHBody &sph_body, std::string variable_name)
			: ParticleDynamicsReduce<VariableType, ReduceSum<VariableType>>(sph_body),
			  GeneralDataDelegateSimple(sph_body),
			  variable_(*particles_->getVariableByName<DataTypeIndex, VariableType>(variable_name))
		{
			this->initial_reference_ = VariableType(0);
		};
		virtual ~BodySummation(){};

	protected:
		StdLargeVec<VariableType> &variable_;
		VariableType ReduceFunction(size_t index_i, Real dt = 0.0) override
		{
			return variable_[index_i];
		};
	};

	/**
	* @class BodyMoment
	* @brief Compute the moment of a body
	*/
	template <int DataTypeIndex, typename VariableType>
	class BodyMoment : public BodySummation<DataTypeIndex, VariableType>
	{
	public:
		explicit BodyMoment(SPHBody &sph_body, std::string variable_name)
			: BodySummation<DataTypeIndex, VariableType>(sph_body, variable_name),
			  mass_(this->particles_->mass_){};
		virtual ~BodyMoment(){};

	protected:
		StdLargeVec<Real> &mass_;
		VariableType ReduceFunction(size_t index_i, Real dt = 0.0) override
		{
			return mass_[index_i] * this->variable_[index_i];
		};
	};

	/**
	 * @class TotalMechanicalEnergy
	 * @brief Compute the total mechanical (kinematic and potential) energy
	 */
	class TotalMechanicalEnergy
		: public ParticleDynamicsReduce<Real, ReduceSum<Real>>,
		  public GeneralDataDelegateSimple
	{
	private:
		UniquePtrKeeper<Gravity> gravity_ptr_keeper_;

	public:
		explicit TotalMechanicalEnergy(SPHBody &sph_body);
		TotalMechanicalEnergy(SPHBody &sph_body, Gravity &gravity_ptr);
		virtual ~TotalMechanicalEnergy(){};

	protected:
		StdLargeVec<Real> &mass_;
		StdLargeVec<Vecd> &vel_n_, &pos_n_;
		Gravity *gravity_;
		Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
	};
}
#endif //GENERAL_DYNAMICS_H