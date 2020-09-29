/**
* @file 	particle_dynamics_bodypart.cpp
* @brief 	This is the implementation of the template class
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#include "particle_dynamics_bodypart.h"

namespace SPH {
	//=================================================================================================//
	void PartDynamicsByParticle::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t i = 0; i < constrained_particles_.size(); ++i)
		{
			Update(constrained_particles_[i], dt);
		}
	}
	//=================================================================================================//
	void PartDynamicsByParticle::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, constrained_particles_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				Update(constrained_particles_[i], dt);
			}
		}, ap);
	}
	//=================================================================================================//
	void PartDynamicsByCell::exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		for (size_t i = 0; i != constrained_cells_.size(); ++i) {
			CellListDataVector& list_data = constrained_cells_[i]->cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num) Update(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void PartDynamicsByCell::parallel_exec(Real dt)
	{
		setBodyUpdated();
		setupDynamics(dt);
		parallel_for(blocked_range<size_t>(0, constrained_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					CellListDataVector& list_data = constrained_cells_[i]->cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num) Update(list_data[num].first, dt);
				}
			}, ap);
	}
	//=================================================================================================//
}