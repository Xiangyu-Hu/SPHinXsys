/**
* @file 	particle_dynamics_constraint.hpp
* @brief 	This is the implementation of the template class
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "particle_dynamics_constraint.h"

namespace SPH {
	//===============================================================//
	template <class BodyType, class ParticlesType, class LagrangianBodyPartType>
	void LagrangianConstraint<BodyType, ParticlesType, LagrangianBodyPartType>::exec(Real dt) 
	{
		PrepareConstraint();
		for (size_t i = 0; i < constrained_particles_.size(); ++i)
		{
			ConstraintAParticle(constrained_particles_[i], dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class LagrangianBodyPartType>
	void LagrangianConstraint<BodyType, ParticlesType, LagrangianBodyPartType>::parallel_exec(Real dt) 
	{
		PrepareConstraint();
		parallel_for(blocked_range<size_t>(0, constrained_particles_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				ConstraintAParticle(constrained_particles_[i], dt);
			}
		}, ap);
	}
	//===============================================================//
	template <class ReturnType, class BodyType, class ParticlesType, class LagrangianBodyPartType>
	ReturnType LagrangianConstraintReduce<ReturnType, BodyType, ParticlesType, LagrangianBodyPartType>
		::exec(Real dt)
	{
		ReturnType temp = initial_reference_;
		SetupReduce();
		//note that base member need to refered by pointer
		//due to the template class has not been instantiated yet
		for (size_t i = 0; i < constrained_particles_.size(); ++i)
		{
			temp = ReduceOperation(temp, ReduceFunction(constrained_particles_[i], dt));
		}
		return OutputResult(temp);
	}
	//===============================================================//
	template <class ReturnType, class BodyType, class ParticlesType, class LagrangianBodyPartType>
	ReturnType LagrangianConstraintReduce<ReturnType, BodyType, ParticlesType, LagrangianBodyPartType>
		::parallel_exec(Real dt)
	{
		ReturnType temp = initial_reference_;
		SetupReduce();
		//note that base member need to refered by pointer
		//due to the template class has not been instantiated yet
		temp = parallel_reduce(blocked_range<size_t>(0, constrained_particles_.size()),
			temp,
			[&](const blocked_range<size_t>& r, ReturnType temp0)->ReturnType {
			for (size_t n = r.begin(); n != r.end(); ++n) {
				temp0 = ReduceOperation(temp0, ReduceFunction(constrained_particles_[n], dt));
			}
			return temp0;
		},
			[this](ReturnType x, ReturnType y)->ReturnType {
			return ReduceOperation(x, y);
		}
		);

		return OutputResult(temp);
	}	
	//===============================================================//
}