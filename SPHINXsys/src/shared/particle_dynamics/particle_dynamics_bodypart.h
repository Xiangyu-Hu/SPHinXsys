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
 * @file 	particle_dynamics_bodypart.h
 * @brief 	Dynamics for bodypart.
 * The dynamics is constrained to a part of the body, 
 * such as in a subregion or on the surface of the body. 
 * The particles of a body part can be defined in an Eulerian or Lagrangian fashion.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_DYNAMICS_BODYPART_H
#define PARTICLE_DYNAMICS_BODYPART_H

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH
{

	/**
	 * @class PartDynamicsByParticle
	 * @brief Abstract class for imposing body part dynamics by particles.
	 * That is the constrained particles will be the same
	 * during the simulation.
	 */
	class PartDynamicsByParticle : public ParticleDynamics<void>
	{
	public:
		PartDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part)
			: ParticleDynamics<void>(sph_body),
			  body_part_particles_(body_part.body_part_particles_){};
		virtual ~PartDynamicsByParticle(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		IndexVector &body_part_particles_;
		ParticleFunctor particle_functor_;
	};

	/**
	 * @class PartSimpleDynamicsByParticle
	 * @brief Abstract class for body part simple particle dynamics.
	 */
	class PartSimpleDynamicsByParticle : public PartDynamicsByParticle
	{
	public:
		PartSimpleDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartSimpleDynamicsByParticle(){};

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
	};

	/**
	 * @class PartInteractionDynamicsByParticle
	 * @brief Abstract class for particle interaction involving in a body part.
	 */
	class PartInteractionDynamicsByParticle : public PartDynamicsByParticle
	{
	public:
		PartInteractionDynamicsByParticle(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartInteractionDynamicsByParticle(){};

	protected:
		virtual void Interaction(size_t index_i, Real dt = 0.0) = 0;
	};

	/**
	 * @class PartInteractionDynamicsByParticleWithUpdate
	 * @brief Abstract class for particle interaction involving in a body part with an extra update step.
	 */
	class PartInteractionDynamicsByParticleWithUpdate : public PartInteractionDynamicsByParticle
	{
	public:
		PartInteractionDynamicsByParticleWithUpdate(SPHBody &sph_body, BodyPartByParticle &body_part);
		virtual ~PartInteractionDynamicsByParticleWithUpdate(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_update_;
	};

	/**
	 * @class PartInteractionDynamicsByParticleWithUpdate
	 * @brief Abstract class for particle interaction involving in a body part with an extra update step.
	 */
	class PartInteractionDynamicsByParticle1Level : public PartInteractionDynamicsByParticleWithUpdate
	{
	public:
		PartInteractionDynamicsByParticle1Level(SPHBody &sph_body, BodyPartByParticle &body_part)
			: PartInteractionDynamicsByParticleWithUpdate(sph_body, body_part),
			  functor_initialization_(std::bind(&PartInteractionDynamicsByParticle1Level::Initialization,
												this, _1, _2)){};
		virtual ~PartInteractionDynamicsByParticle1Level(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		virtual void Initialization(size_t index_i, Real dt = 0.0) = 0;
		ParticleFunctor functor_initialization_;
	};

	/**
	 * @class PartDynamicsByCell
	 * @brief Abstract class for imposing Eulerian constrain to a body.
	 * The constrained particles are in the tagged cells .
	 */
	class PartDynamicsByCell : public ParticleDynamics<void>
	{
	public:
		PartDynamicsByCell(SPHBody &sph_body, BodyPartByCell &body_part)
			: ParticleDynamics<void>(sph_body),
			  body_part_cells_(body_part.body_part_cells_){};
		virtual ~PartDynamicsByCell(){};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;

	protected:
		CellLists &body_part_cells_;

		virtual void Update(size_t index_i, Real dt = 0.0) = 0;
	};
	/**
	  * @class PartDynamicsByCellReduce
	  * @brief Abstract class for reduce operation in a Eulerian constrain region.
	  */
	template <class ReturnType, typename ReduceOperation>
	class PartDynamicsByCellReduce : public ParticleDynamics<ReturnType>
	{
	public:
		PartDynamicsByCellReduce(SPHBody &sph_body, BodyPartByCell &body_part)
			: ParticleDynamics<ReturnType>(sph_body), body_part_cells_(body_part.body_part_cells_),
			  quantity_name_("ReducedQuantity"), initial_reference_(){};
		virtual ~PartDynamicsByCellReduce(){};

		ReturnType InitialReference() { return initial_reference_; };
		std::string QuantityName() { return quantity_name_; };

		virtual ReturnType exec(Real dt = 0.0) override
		{
			ReturnType temp = initial_reference_;
			this->SetupReduce();
			for (size_t i = 0; i != body_part_cells_.size(); ++i)
			{
				ListDataVector &list_data = body_part_cells_[i]->cell_list_data_;
				for (size_t num = 0; num < list_data.size(); ++num)
				{
					temp = reduce_operation_(temp, ReduceFunction(list_data[num].first, dt));
				}
			}
			return OutputResult(temp);
		};

		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			ReturnType temp = initial_reference_;
			this->SetupReduce();
			temp = parallel_reduce(
				blocked_range<size_t>(0, body_part_cells_.size()),
				temp,
				[&](const blocked_range<size_t> &r, ReturnType temp0) -> ReturnType
				{
					for (size_t i = r.begin(); i != r.end(); ++i)
					{
						ListDataVector &list_data = body_part_cells_[i]->cell_list_data_;
						for (size_t num = 0; num < list_data.size(); ++num)
						{
							temp0 = reduce_operation_(temp0, ReduceFunction(list_data[num].first, dt));
						}
					}
					return temp0;
				},
				[this](ReturnType x, ReturnType y) -> ReturnType
				{ return reduce_operation_(x, y); });

			return OutputResult(temp);
		};

	protected:
		ReduceOperation reduce_operation_;
		CellLists &body_part_cells_;
		std::string quantity_name_;
		ReturnType initial_reference_;
		virtual void SetupReduce(){};
		virtual ReturnType ReduceFunction(size_t index_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
	};
	/** 
	  * @class PartDynamicsByParticleReduce
	  * @brief reduce operation in a Lagrangian contrained region.
	  */
	template <class ReturnType, typename ReduceOperation>
	class PartDynamicsByParticleReduce : public ParticleDynamics<ReturnType>
	{
	public:
		PartDynamicsByParticleReduce(SPHBody &sph_body, BodyPartByParticle &body_part)
			: ParticleDynamics<ReturnType>(sph_body),
			  body_part_particles_(body_part.body_part_particles_),
			  quantity_name_("ReducedQuantity"), initial_reference_(){};
		virtual ~PartDynamicsByParticleReduce(){};

		ReturnType InitialReference() { return initial_reference_; };
		std::string QuantityName() { return quantity_name_; };

		virtual ReturnType exec(Real dt = 0.0) override
		{
			ReturnType temp = initial_reference_;
			this->SetupReduce();
			for (size_t i = 0; i < body_part_particles_.size(); ++i)
			{
				temp = reduce_operation_(temp, ReduceFunction(body_part_particles_[i], dt));
			}
			return OutputResult(temp);
		};
		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			ReturnType temp = initial_reference_;
			this->SetupReduce();
			temp = parallel_reduce(
				blocked_range<size_t>(0, body_part_particles_.size()),
				temp,
				[&](const blocked_range<size_t> &r, ReturnType temp0) -> ReturnType
				{
					for (size_t n = r.begin(); n != r.end(); ++n)
					{
						temp0 = reduce_operation_(temp0, ReduceFunction(body_part_particles_[n], dt));
					}
					return temp0;
				},
				[this](ReturnType x, ReturnType y) -> ReturnType
				{
					return reduce_operation_(x, y);
				});

			return OutputResult(temp);
		};

	protected:
		ReduceOperation reduce_operation_;
		IndexVector &body_part_particles_;
		std::string quantity_name_;
		ReturnType initial_reference_;
		virtual void SetupReduce(){};
		virtual ReturnType ReduceFunction(size_t index_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
	};
}
#endif //PARTICLE_DYNAMICS_BODYPART_H