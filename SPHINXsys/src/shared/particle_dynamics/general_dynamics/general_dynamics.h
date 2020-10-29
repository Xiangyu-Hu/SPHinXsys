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
	* set the number of ghost particles into zero.
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
	 * @class PeriodicConditionInAxisDirection
	 * @brief Base class for two different type periodic boundary conditions.
	 */
	class PeriodicConditionInAxisDirection : public BoundingInAxisDirection
	{
	protected:
		/** cells in which particle checked for bounding */
		StdVec<CellVector> bound_cells_;

		/**
		* @class PeriodicBounding
		* @brief Periodic bounding particle position in an axis direction
		*/
		class PeriodicBounding : public BoundingInAxisDirection
		{
		protected:
			Vecd periodic_translation_;
			//cells in which particle checked for bounding
			StdVec<CellVector>& bound_cells_;

			/**compute the distance for periodic translation. */
			void setPeriodicTranslation();

			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);
		public:
			PeriodicBounding(StdVec<CellVector>& bound_cells, SPHBody* body, int axis_direction);
			virtual ~PeriodicBounding() {};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};
	public:
		PeriodicConditionInAxisDirection(SPHBody* body, int axis_direction);
		virtual ~PeriodicConditionInAxisDirection() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/**
	 * @class PeriodicConditionInAxisDirectionUsingCellLinkedList
	 * @brief The method imposing periodic boundary condition in an axis direction.
	 *	It includes two different steps, i.e. imposing periodic bounding and condition.
	 *	The first step is carried out before update cell linked list and
	 *	the second after the updating.
	 *	If the exec or parallel_exec is called directly, error message will be given.
	 */
	class PeriodicConditionInAxisDirectionUsingCellLinkedList : 
		public PeriodicConditionInAxisDirection
	{
	protected:
		/**
		* @class PeriodicCondition
		* @brief Periodic boundary condition in an axis direction
		*/
		class PeriodicCellLinkedList : public PeriodicBounding
		{
		protected:
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;
		public:

			PeriodicCellLinkedList(StdVec<CellVector>& bound_cells, SPHBody* body, int axis_direction)
				: PeriodicBounding(bound_cells, body, axis_direction) {};
			virtual ~PeriodicCellLinkedList() {};

			/** This class is only implemented in sequential due to memory conflicts.
			 * Because the cell list data is not concurrent vector.
			 */
			virtual void parallel_exec(Real dt = 0.0) override { exec(); };
		};
	public:
		PeriodicConditionInAxisDirectionUsingCellLinkedList(SPHBody* body, int axis_direction) :
			PeriodicConditionInAxisDirection(body, axis_direction),
			bounding_(this->bound_cells_, body, axis_direction),
			update_cell_linked_list_(this->bound_cells_, body, axis_direction) {};
		virtual ~PeriodicConditionInAxisDirectionUsingCellLinkedList() {};

		PeriodicBounding bounding_;
		PeriodicCellLinkedList update_cell_linked_list_;
	};

	/**
	 * @class PeriodicConditionInAxisDirectionUsingGhostParticles
	 * @brief The method imposing periodic boundary condition in an axis direction by using ghost particles.
	 *	It includes three different steps, i.e. imposing periodic bounding, creating ghosts and update ghost state.
	 *	The first step is carried out before update cell linked list and
	 *	the second and third after the updating.
	 *	If the exec or parallel_exec is called directly, error message will be given.
	 */
	class PeriodicConditionInAxisDirectionUsingGhostParticles : 
		public PeriodicConditionInAxisDirection
	{
	protected:
		/** ghost particles createded for impose boundary condition. */
		StdVec<IndexVector> ghost_particles_;

		/**
		 * @class CreatPeriodicGhostParticles
		 * @brief create ghost particles in an axis direction
		 */
		class CreatPeriodicGhostParticles : public PeriodicBounding
		{
		protected:
			StdVec<IndexVector>& ghost_particles_;
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0) override;
		public:
			CreatPeriodicGhostParticles(StdVec<CellVector>& bound_cells,
				StdVec<IndexVector>& ghost_particles, SPHBody* body, int axis_direction) :
				PeriodicBounding(bound_cells, body, axis_direction),
				ghost_particles_(ghost_particles) {};
			virtual ~CreatPeriodicGhostParticles() {};

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
			StdVec<IndexVector>& ghost_particles_;
			void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			void checkUpperBound(size_t index_i, Real dt = 0.0) override;
		public:
			UpdatePeriodicGhostParticles(StdVec<CellVector>& bound_cells,
				StdVec<IndexVector>& ghost_particles, SPHBody* body, int axis_direction) :
				PeriodicBounding(bound_cells, body, axis_direction), ghost_particles_(ghost_particles) {};
			virtual ~UpdatePeriodicGhostParticles() {};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};
	public:
		PeriodicConditionInAxisDirectionUsingGhostParticles(SPHBody* body, int axis_direction) :
			PeriodicConditionInAxisDirection(body, axis_direction),
			bounding_(this->bound_cells_, body, axis_direction),
			ghost_creation_(this->bound_cells_, this->ghost_particles_, body, axis_direction),
			ghost_update_(this->bound_cells_, this->ghost_particles_, body, axis_direction)
		{
			ghost_particles_.resize(2);
		};

		virtual ~PeriodicConditionInAxisDirectionUsingGhostParticles() {};

		PeriodicBounding bounding_;
		CreatPeriodicGhostParticles ghost_creation_;
		UpdatePeriodicGhostParticles ghost_update_;
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

		class MirrorBounding : public BoundingInAxisDirection
		{
		protected:
			CellVector& bound_cells_;
			virtual void checkLowerBound(size_t index_i, Real dt = 0.0);
			virtual void checkUpperBound(size_t index_i, Real dt = 0.0);
			ParticleFunctor checking_bound_;

			StdLargeVec<Vecd>& vel_n_;
			/** mirror the particle physical state along an axis direction. */
			void mirrorInAxisDirection(size_t particle_index_i, Vecd body_bound, int axis_direction);
		public:
			MirrorBounding(CellVector& bound_cells, SPHBody* body, int axis_direction, bool positive);
			virtual ~MirrorBounding() {};
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
		class UpdatingGhostStates : public MirrorBounding
		{
		protected:
			IndexVector& ghost_particles_;
			void checkLowerBound(size_t index_i, Real dt = 0.0) override;
			void checkUpperBound(size_t index_i, Real dt = 0.0) override;
			ParticleFunctor checking_bound_update_;
		public:
			UpdatingGhostStates(IndexVector& ghost_particles, CellVector& bound_cells,
				SPHBody* body, int axis_direction, bool positive);
			virtual ~UpdatingGhostStates() {};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

	public:
		MirrorBounding bounding_;
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

	/**
	 * @class DampingBySplittingAlgorithm
	 * @brief A quantity damping by splitting scheme
	 * this method modifies the quantity directly.
	 * Note that, if periodic boundary condition is applied, 
	 * the parallelized version of the method requires the one using ghost particles 
	 * because the splitting partition only works in this case.  
	 */
	template<typename VariableType>
	class DampingBySplittingAlgorithm : 
		public InteractionDynamicsSplitting,
		public DataDelegateInner<SPHBody, BaseParticles, BaseMaterial>
	{
	public:
		DampingBySplittingAlgorithm(SPHBodyInnerRelation* body_inner_relation,
			StdLargeVec<VariableType>& variable, Real eta) :
			InteractionDynamicsSplitting(body_inner_relation->sph_body_),
			DataDelegateInner<SPHBody, BaseParticles, BaseMaterial>(body_inner_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_), variable_(variable), eta_(eta) {};
		virtual ~DampingBySplittingAlgorithm() {};
	protected:
		StdLargeVec<Real>& Vol_, & mass_;
		StdLargeVec<VariableType>& variable_;
		Real eta_; /**< damping coefficient */
		virtual void Interaction(size_t index_i, Real dt = 0.0) override
		{
			Real Vol_i = Vol_[index_i];
			Real mass_i = mass_[index_i];
			VariableType& variable_i = variable_[index_i];

			VariableType error(0);
			Real parameter_a(0);
			Real parameter_c(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				//linear projection 
				VariableType variable_derivative = (variable_i - variable_[index_j]);
				Real parameter_b = 2.0 * eta_ * inner_neighborhood.dW_ij_[n]
					* Vol_i * Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

				error -= variable_derivative * parameter_b;
				parameter_a += parameter_b;
				parameter_c += parameter_b * parameter_b;
			}

			parameter_a -= mass_i;
			Real parameter_l = parameter_a * parameter_a + parameter_c;
			VariableType parameter_k = error / (parameter_l + TinyReal);
			variable_[index_i] += parameter_k * parameter_a;

			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				Real parameter_b = 2.0 * eta_ * inner_neighborhood.dW_ij_[n]
					* Vol_i * Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

				//predicted quantity at particle j
				VariableType variable_j = variable_[index_j] - parameter_k * parameter_b;
				VariableType variable_derivative = (variable_i - variable_j);

				//exchange in conservation form
				variable_[index_j] -= variable_derivative * parameter_b / mass_[index_j];
			}
		};
	};

	/**
	* @class DampingBySplittingPairwise
	* @brief A quantity damping by a pairwise splitting scheme
	* this method modifies the quantity directly
	* Note that, if periodic boundary condition is applied, 
	* the parallelized version of the method requires the one using ghost particles 
	* because the splitting partition only works in this case.  
	*/
	template<typename VariableType>
	class DampingBySplittingPairwise : 
		public InteractionDynamicsSplitting,
		public DataDelegateInner<SPHBody, BaseParticles, BaseMaterial>
	{
	public:
		DampingBySplittingPairwise(SPHBodyInnerRelation* body_inner_relation,
			StdLargeVec<VariableType>& variable, Real eta) :
			InteractionDynamicsSplitting(body_inner_relation->sph_body_),
			DataDelegateInner<SPHBody, BaseParticles, BaseMaterial>(body_inner_relation),
			Vol_(particles_->Vol_), mass_(particles_->mass_), variable_(variable), eta_(eta) {};
		virtual ~DampingBySplittingPairwise() {};
	protected:
		StdLargeVec<Real>& Vol_, & mass_;
		StdLargeVec<VariableType>& variable_;
		Real eta_; /**< damping coefficient */
		
		virtual void Interaction(size_t index_i, Real dt = 0.0) override 
		{
			Real Vol_i = Vol_[index_i];
			Real mass_i = mass_[index_i];
			VariableType& variable_i = variable_[index_i];

			StdVec<Real> parameter_b(MaximumNeighborhoodSize);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			//forward sweep
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real mass_j = mass_[index_j];

				VariableType variable_derivative = (variable_i - variable_[index_j]);
				parameter_b[n] = eta_ * inner_neighborhood.dW_ij_[n]
					* Vol_i * Vol_[index_j] * dt / inner_neighborhood.r_ij_[n];

				VariableType increment = parameter_b[n] * variable_derivative
					/ (mass_i * mass_j - parameter_b[n] * (mass_i + mass_j));
				variable_[index_i] += increment * mass_j;
				variable_[index_j] -= increment * mass_i;
			}

			//backward sweep
			for (size_t n = inner_neighborhood.current_size_; n != 0; --n)
			{
				size_t index_j = inner_neighborhood.j_[n - 1];
				Real mass_j = mass_[index_j];

				VariableType variable_derivative = (variable_i - variable_[index_j]);
				VariableType increment = parameter_b[n - 1] * variable_derivative
					/ (mass_i * mass_j - parameter_b[n - 1] * (mass_i + mass_j));

				variable_[index_i] += increment * mass_j;
				variable_[index_j] -= increment * mass_i;
			}
		};
	};

	/**
	* @class DampingBySplittingWithRandomChoice
	* @brief A random choice method for obstaining static equilibrium state
	* Note that, if periodic boundary condition is applied, 
	* the parallelized version of the method requires the one using ghost particles 
	* because the splitting partition only works in this case.  
	*/
	template<class DampingAlgorithmType, typename VariableType>
	class DampingBySplittingWithRandomChoice : public DampingAlgorithmType
	{
	protected:
		Real random_ratio_;
		bool RandomChoice() 
		{
			return ((double)rand() / (RAND_MAX)) < random_ratio_ ? true : false;
		};
	public:
		DampingBySplittingWithRandomChoice(SPHBodyInnerRelation* body_inner_relation,
			Real random_ratio, StdLargeVec<VariableType>& variable, Real eta) :
			DampingAlgorithmType(body_inner_relation, variable, eta / random_ratio),
			random_ratio_(random_ratio) {};
		virtual ~DampingBySplittingWithRandomChoice() {};

		virtual void exec(Real dt = 0.0) override
		{
			if (RandomChoice()) DampingAlgorithmType::exec(dt);
		};
		virtual void parallel_exec(Real dt = 0.0) override
		{
			if (RandomChoice()) DampingAlgorithmType::parallel_exec(dt);
		};
	};
}
