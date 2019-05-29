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
	template <class BodyType, class ParticlesType, class LagrangianBodyPartType>
	class LagrangianConstraint : public ParticleDynamics<void, BodyType, ParticlesType>
	{
	protected:
		LagrangianBodyPartType *body_part_;
		IndexVector &constrained_particles_;

		virtual void  PrepareConstraint() {};
		virtual void  ConstraintAParticle(size_t index_particle_i, 	Real dt = 0.0) = 0;
	public:
		LagrangianConstraint(BodyType *body, LagrangianBodyPartType *body_part)
			: ParticleDynamics<void, BodyType, ParticlesType>(body),
			body_part_(body_part), constrained_particles_(body_part->body_part_particles_) {};
		virtual ~LagrangianConstraint() {};

		virtual void exec(Real dt = 0.0);
		virtual void parallel_exec(Real dt = 0.0);
	};

	/** 
	 * @class EulerianConstraint
	 * @brief Imposing Eulerian constrain to a body.
	 * The constrained particles are the cells tagged.
	 */
	template <class BodyType, class ParticlesType, class EulerianBodyPartType>
	class EulerianConstraint : public ParticleDynamicsByCells<BodyType, ParticlesType>
	{
	protected:
		EulerianBodyPartType *body_part_;
		CellVector &constrained_cells_;

		virtual void PrepareConstraint() {};
		virtual void  ConstraintAParticle(size_t index_particle_i,
			Real dt = 0.0) = 0;
	public:
		EulerianConstraint(BodyType *body, EulerianBodyPartType *body_part)
			: ParticleDynamicsByCells<BodyType, ParticlesType>(body),
			body_part_(body_part), constrained_cells_(body_part->body_part_cells_) {};
		virtual ~EulerianConstraint() {};

		virtual void exec(Real dt = 0.0) override;
		virtual void parallel_exec(Real dt = 0.0) override;
	};
	/** 
	  * @class LagrangianConstraintReduce
	  * @brief reduce operation in a Lagrangian contrained region.
	  */
	template <class ReturnType, class BodyType, class ParticlesType, class LagrangianBodyPartType>
	class LagrangianConstraintReduce
		: public ParticleDynamicsReduce<ReturnType, BodyType, ParticlesType>
	{
	protected:
		LagrangianBodyPartType *body_part_;
		IndexVector &constrained_particles_;

	public:
		LagrangianConstraintReduce(BodyType* body, LagrangianBodyPartType *body_part)
			: ParticleDynamicsReduce<ReturnType, BodyType, ParticlesType>(body),
			body_part_(body_part), constrained_particles_(body_part->body_part_particles_) {};
		virtual ~LagrangianConstraintReduce() {};

		/** sequential */
		virtual ReturnType exec(Real dt = 0.0);

		/** parallel */
		virtual ReturnType parallel_exec(Real dt = 0.0);
	};

	/** 
	 * @class LagrangianConstraintSum
	 * @brief compute reduced sum form Lagrangian contrained region.
	 */
	template <class ReturnType, class BodyType, class ParticlesType, class LagrangianBodyPartType>
	class LagrangianConstraintSum
		: public LagrangianConstraintReduce<ReturnType, BodyType, ParticlesType, LagrangianBodyPartType>
	{
	protected:

		virtual ReturnType ReduceOperation(ReturnType x, ReturnType y) override { return x + y; };
	public:
		LagrangianConstraintSum(BodyType* body, LagrangianBodyPartType *body_part)
			: LagrangianConstraintReduce<ReturnType, BodyType, ParticlesType, LagrangianBodyPartType>(body, body_part) {};
		virtual ~LagrangianConstraintSum() {};
	};
}
