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
 * @file 	particle_dynamics_constraint.h
 * @brief 	This is the class for constrain bodies.
 * We constrain the particles on the body. These particles can be
 * in a subregion of on the surface of the body. 
 * We also consider Eulerian and Lagrangian constraint approaches
 * @author	Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_particle_dynamics.h"
#include "base_particle_dynamics.hpp"

namespace SPH {

	/**
	 * @class PartDynamicsByParticle
	 * @brief Imposing Lagrangian constrain to a body.
	 * That is the constrained particles will be the same
	 * during the simulation.
	 */
	class PartDynamicsByParticle : public ParticleDynamics<void>
	{
	public:
		PartDynamicsByParticle(SPHBody* sph_body, BodyPartByParticle *body_part)
			: ParticleDynamics<void>(sph_body),
			constrained_particles_(body_part->body_part_particles_) {};
		virtual ~PartDynamicsByParticle() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	protected:
		IndexVector& constrained_particles_;

		virtual void  Update(size_t index_i, Real dt = 0.0) = 0;
	};

	/** 
	 * @class PartDynamicsByCell
	 * @brief Imposing Eulerian constrain to a body.
	 * The constrained particles are in the tagged cells .
	 */
	class PartDynamicsByCell : public ParticleDynamics<void>
	{
	public:
		PartDynamicsByCell(SPHBody* sph_body, BodyPartByCell *body_part)
			: ParticleDynamics<void>(sph_body),
			constrained_cells_(body_part->body_part_cells_) {};
		virtual ~PartDynamicsByCell() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	protected:
		CellLists& constrained_cells_;

		virtual void  Update(size_t index_i, Real dt = 0.0) = 0;
	};
	/** 
	  * @class PartDynamicsByParticleReduce
	  * @brief reduce operation in a Lagrangian contrained region.
	  */
	template <class ReturnType, typename ReduceOperation>
	class PartDynamicsByParticleReduce : public ParticleDynamics<ReturnType>
	{
	public:
		PartDynamicsByParticleReduce(SPHBody* sph_body, BodyPartByParticle *body_part)
			: ParticleDynamics<ReturnType>(sph_body),
			constrained_particles_(body_part->body_part_particles_) {};
		virtual ~PartDynamicsByParticleReduce() {};

		virtual ReturnType exec(Real dt = 0.0) override
		{
			ReturnType temp = initial_reference_;
			this->SetupReduce();
			//note that base member need to referred by pointer
			//due to the template class has not been instantiated yet
			for (size_t i = 0; i < constrained_particles_.size(); ++i)
			{
				temp = reduce_operation_(temp, ReduceFunction(constrained_particles_[i], dt));
			}
			return OutputResult(temp);
		};
		virtual ReturnType parallel_exec(Real dt = 0.0) override
		{
			ReturnType temp = initial_reference_;
			this->SetupReduce();
			//note that base member need to referred by pointer
			//due to the template class has not been instantiated yet
			temp = parallel_reduce(blocked_range<size_t>(0, constrained_particles_.size()),
				temp,
				[&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType {
					for (size_t n = r.begin(); n != r.end(); ++n) {
						temp0 = reduce_operation_(temp0, ReduceFunction(constrained_particles_[n], dt));
					}
					return temp0;
				},
				[this](ReturnType x, ReturnType y)->ReturnType {
					return reduce_operation_(x, y);
				}
				);

			return OutputResult(temp);
		};
	protected:
		ReduceOperation reduce_operation_;

		IndexVector& constrained_particles_;

		//inital or reference value
		ReturnType initial_reference_;
		virtual void SetupReduce() {};
		virtual ReturnType ReduceFunction(size_t index_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
	};
}
