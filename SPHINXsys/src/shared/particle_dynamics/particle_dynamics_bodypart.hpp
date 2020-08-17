/**
* @file 	particle_dynamics_constraint.hpp
* @brief 	This is the implementation of the template class
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "particle_dynamics_bodypart.h"

namespace SPH {
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType>
	void PartDynamicsByParticle<BodyType, ParticlesType, BodyPartByParticleType, MaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		for (size_t i = 0; i < constrained_particles_.size(); ++i)
		{
			Update(constrained_particles_[i], dt);
		}
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType>
	void PartDynamicsByParticle<BodyType, ParticlesType, BodyPartByParticleType, MaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, constrained_particles_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Update(constrained_particles_[i], dt);
			}
		}, ap);
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class BodyPartByCellType, class MaterialType>
	void PartDynamicsByCell<BodyType, ParticlesType, BodyPartByCellType, MaterialType>
		::exec(Real dt)
	{
		this->setupDynamics(dt);
		for (size_t i = 0; i != constrained_cells_.size(); ++i) {
			CellListDataVector& list_data = constrained_cells_[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num) Update(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	template <class BodyType, class ParticlesType, class BodyPartByCellType, class MaterialType>
	void PartDynamicsByCell<BodyType, ParticlesType, BodyPartByCellType, MaterialType>
		::parallel_exec(Real dt)
	{
		this->setupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, constrained_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					CellListDataVector& list_data = constrained_cells_[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num) Update(list_data[num].first, dt);
				}
			}, ap);
	}
	//=================================================================================================//
	template <class ReturnType, typename ReduceOperation,
		class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType>
		ReturnType PartDynamicsByParticleReduce<ReturnType, ReduceOperation,
		BodyType, ParticlesType, BodyPartByParticleType, MaterialType>
		::exec(Real dt)
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
	}
	//=================================================================================================//
	template <class ReturnType, class ReduceOperation,
		class BodyType, class ParticlesType, class BodyPartByParticleType, class MaterialType>
		ReturnType PartDynamicsByParticleReduce<ReturnType, ReduceOperation,
		BodyType, ParticlesType, BodyPartByParticleType, MaterialType>
		::parallel_exec(Real dt)
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
	}	
	//=================================================================================================//
}