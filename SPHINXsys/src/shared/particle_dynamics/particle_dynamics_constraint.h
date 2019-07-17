/**
 * @file 	particle_dynamics_constraint.h
 * @brief 	This is the class for constrain bodies.
 * We constrain the particles on the body. These paericles can be
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
	 * @class LagrangianConstraint
	 * @brief Imposing Lagrangian constrain to a body.
	 * That is the constrained particles will be the same
	 * during the simulation.
	 */
	template <class BodyType, class ParticlesType, class LagrangianBodyPartType, class MaterialType = Material>
	class LagrangianConstraint : public ParticleDynamics<void, BodyType, ParticlesType, MaterialType>
	{
	protected:
		LagrangianBodyPartType *body_part_;
		IndexVector &constrained_particles_;

		virtual void  PrepareConstraint() {};
		virtual void  ConstraintAParticle(size_t index_particle_i, 	Real dt = 0.0) = 0;
	public:
		LagrangianConstraint(BodyType *body, LagrangianBodyPartType *body_part)
			: ParticleDynamics<void, BodyType, ParticlesType, MaterialType>(body),	body_part_(body_part), 
			constrained_particles_(body_part->body_part_particles_) {};
		virtual ~LagrangianConstraint() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};

	/** 
	 * @class EulerianConstraint
	 * @brief Imposing Eulerian constrain to a body.
	 * The constrained particles are the cells tagged.
	 */
	template <class BodyType, class ParticlesType, class EulerianBodyPartType, class MaterialType = Material>
	class EulerianConstraint : public ParticleDynamicsByCells<BodyType, ParticlesType, MaterialType >
	{
	protected:
		EulerianBodyPartType *body_part_;
		CellVector &constrained_cells_;

		virtual void PrepareConstraint() {};
		virtual void  ConstraintAParticle(size_t index_particle_i,
			Real dt = 0.0) = 0;
	public:
		EulerianConstraint(BodyType *body, EulerianBodyPartType *body_part)
			: ParticleDynamicsByCells<BodyType, ParticlesType, MaterialType>(body), 
			body_part_(body_part), constrained_cells_(body_part->body_part_cells_) {};
		virtual ~EulerianConstraint() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
	/** 
	  * @class LagrangianConstraintReduce
	  * @brief reduce operation in a Lagrangian contrained region.
	  */
	template <class ReturnType, typename ReduceOperation, 
		class BodyType, class ParticlesType, class LagrangianBodyPartType, class MaterialType = Material>
	class LagrangianConstraintReduce : public ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>
	{
	protected:
		ReduceOperation reduce_operation_;

		LagrangianBodyPartType *body_part_;
		IndexVector &constrained_particles_;

		//inital or refence value
		ReturnType initial_reference_;
		virtual void SetupReduce() {};
		virtual ReturnType ReduceFunction(size_t index_particle_i, Real dt = 0.0) = 0;
		virtual ReturnType OutputResult(ReturnType reduced_value) { return reduced_value; };
	public:
		LagrangianConstraintReduce(BodyType* body, LagrangianBodyPartType *body_part)
			: ParticleDynamics<ReturnType, BodyType, ParticlesType, MaterialType>(body), body_part_(body_part), 
			constrained_particles_(body_part->body_part_particles_) {};
		virtual ~LagrangianConstraintReduce() {};

		/** sequential */
		virtual ReturnType exec(Real dt = 0.0) override;

		/** parallel */
		virtual ReturnType parallel_exec(Real dt = 0.0) override;
	};
}
