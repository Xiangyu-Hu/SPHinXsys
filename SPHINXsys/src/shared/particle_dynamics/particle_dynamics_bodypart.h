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
	template <class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType = BaseMaterial>
	class PartDynamicsByParticle : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		BodyPartByParticleType *body_part_;
		IndexVector &constrained_particles_;

		virtual void  Update(size_t index_particle_i, 	Real dt = 0.0) = 0;
	public:
		PartDynamicsByParticle(SPHBody* body, BodyPartByParticleType *body_part)
			: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body),	body_part_(body_part), 
			constrained_particles_(body_part->body_part_particles_) {};
		virtual ~PartDynamicsByParticle() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/** 
	 * @class PartDynamicsByCell
	 * @brief Imposing Eulerian constrain to a body.
	 * The constrained particles are the cells tagged.
	 */
	template <class BodyType, class ParticlesType, class BodyPartByCellType, class MaterialType = BaseMaterial>
	class PartDynamicsByCell : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType >
	{
	protected:
		BodyPartByCellType *body_part_;
		CellLists &constrained_cells_;

		virtual void  Update(size_t index_particle_i, Real dt = 0.0) = 0;
	public:
		PartDynamicsByCell(SPHBody* body, BodyPartByCellType *body_part)
			: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body), 
			body_part_(body_part), constrained_cells_(body_part->body_part_cells_) {};
		virtual ~PartDynamicsByCell() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
	/** 
	  * @class PartDynamicsByParticleReduce
	  * @brief reduce operation in a Lagrangian contrained region.
	  */
	template <class ReturnType, typename ReduceOperation, 
		class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType = BaseMaterial>
	class PartDynamicsByParticleReduce : public ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>
	{
	protected:
		ReduceOperation reduce_operation_;

		BodyPartByParticleType *body_part_;
		IndexVector &constrained_particles_;

		//inital or reference value
		ReturnType initial_reference_;
		virtual void SetupReduce() {};
		virtual ReturnType ReduceFunction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
	public:
		PartDynamicsByParticleReduce(BodyType* body, BodyPartByParticleType *body_part)
			: ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>(body), body_part_(body_part), 
			constrained_particles_(body_part->body_part_particles_) {};
		virtual ~PartDynamicsByParticleReduce() {};

		virtual ReturnType exec(Real dt = 0.0) override;
		virtual ReturnType parallel_exec(Real dt = 0.0) override;
	};
}
