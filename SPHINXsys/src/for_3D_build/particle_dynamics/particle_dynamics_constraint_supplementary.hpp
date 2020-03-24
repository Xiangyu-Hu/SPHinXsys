/**
* @file 	particle_dynamics_constraint.hpp
* @brief 	This is the implementation of the template class for 3D build
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "particle_dynamics_constraint.h"

namespace SPH {
	//===============================================================//
	template <class BodyType, class ParticlesType, class BodyPartByCellType, class MaterialType>
	void ConstraintByCell<BodyType, ParticlesType, BodyPartByCellType, MaterialType>
		::exec(Real dt)
	{
		PrepareConstraint();

		for (size_t i = 0; i != constrained_cells_.size(); ++i) {
			ConcurrentListDataVector &list_data
				= this->cell_linked_lists_[constrained_cells_[i][0]][constrained_cells_[i][1]][constrained_cells_[i][2]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				ConstraintAParticle(list_data[num].first, dt);
		}
	}
	//===============================================================//
	template <class BodyType, class ParticlesType, class BodyPartByCellType, class MaterialType>
	void ConstraintByCell<BodyType, ParticlesType, BodyPartByCellType, MaterialType>
		::parallel_exec(Real dt)
	{
		PrepareConstraint();

		parallel_for(blocked_range<size_t>(0, constrained_cells_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				ConcurrentListDataVector &list_data
					= this->cell_linked_lists_[constrained_cells_[i][0]][constrained_cells_[i][1]][constrained_cells_[i][2]]
					.particle_data_lists_;
				for (size_t num = 0; num < list_data.size(); ++num)
					ConstraintAParticle(list_data[num].first, dt);
			}
		}, ap);
	}
	//===============================================================//
}