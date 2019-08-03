/**
 * @file 	general_dynamics_supplementary.cpp
 * @brief 	supplementary functions declared in general_dynamics.h.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			this code developed based one the beta version
 * 			add Bounding particles in y direction;
 * 			add Periodic condition in y directin;
 * 			-- Chi
 */
#include "general_dynamics.h"

namespace SPH 
{
//=================================================================================================//
	BoundingInXDirection::BoundingInXDirection(SPHBody* body)
		: BoundingBodyDomain(body)
	{
		//lower bound cells
		for (size_t j = SMAX(int(body_lower_bound_cell_[1] - 1), 0);
				j < SMIN(int(body_upper_bound_cell_[1] + 2), int(number_of_cells_[1])); ++j)
			for (size_t i = SMAX(int(body_lower_bound_cell_[0]) - 1, 0);
					i <= SMIN(int(body_lower_bound_cell_[0] + 1), int(number_of_cells_[0] - 1)); ++i) {
			lower_bound_cells_.push_back(Vecu(i, j));
		}

		//upper bound cells
		for (size_t j = SMAX(int(body_lower_bound_cell_[1] - 1), 0);
			j < SMIN(int(body_upper_bound_cell_[1] + 2), int(number_of_cells_[1])); ++j)
			for (size_t i = SMAX(int(body_upper_bound_cell_[0]) - 1, 0);
				i <= SMIN(int(body_upper_bound_cell_[0] + 1), int(number_of_cells_[0] - 1)); ++i) {
			upper_bound_cells_.push_back(Vecu(i, j));
		}
	}
//=================================================================================================//
	void BoundingInXDirection::exec(Real dt)
	{
		//check lower bound
		for (size_t i = 0; i != lower_bound_cells_.size(); ++i) {
			ListDataVector &list_data
				= cell_linked_lists_[lower_bound_cells_[i][0]][lower_bound_cells_[i][1]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				CheckLowerBound(list_data[num].first, list_data[num].second, dt);
		}

		//check upper bound
		for (size_t i = 0; i != upper_bound_cells_.size(); ++i) {
			ListDataVector &list_data
				= cell_linked_lists_[upper_bound_cells_[i][0]][upper_bound_cells_[i][1]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				CheckUpperBound(list_data[num].first, list_data[num].second, dt);
		}
	}
//=================================================================================================//
	void BoundingInXDirection::parallel_exec(Real dt)
	{
		//check lower bound
		parallel_for(blocked_range<size_t>(0, lower_bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				ListDataVector &list_data
					= cell_linked_lists_[lower_bound_cells_[i][0]][lower_bound_cells_[i][1]]
					.particle_data_lists_;
				for (size_t num = 0; num < list_data.size(); ++num)
					CheckLowerBound(list_data[num].first, list_data[num].second, dt);
			}
		}, ap);

		//check upper bound
		parallel_for(blocked_range<size_t>(0, upper_bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				ListDataVector &list_data
					= cell_linked_lists_[upper_bound_cells_[i][0]][upper_bound_cells_[i][1]]
					.particle_data_lists_;
				for (size_t num = 0; num < list_data.size(); ++num)
					CheckUpperBound(list_data[num].first, list_data[num].second, dt);
			}
		}, ap);
	}
//=================================================================================================//
	BoundingInYDirection::BoundingInYDirection(SPHBody* body)
		: BoundingBodyDomain(body)
	{
		//lower bound cells
		for (size_t i = SMAX(int(body_lower_bound_cell_[0] - 1), 0);
				i < SMIN(int(body_upper_bound_cell_[0] + 2), int(number_of_cells_[0])); ++i)
			for (size_t j = SMAX(int(body_lower_bound_cell_[1]) - 1, 0);
					j <= SMIN(int(body_lower_bound_cell_[1] + 1), int(number_of_cells_[1] - 1)); ++j) 
			{
				lower_bound_cells_.push_back(Vecu(i, j));
			}

		//upper bound cells
		for (size_t i = SMAX(int(body_lower_bound_cell_[0] - 1), 0);
				i < SMIN(int(body_upper_bound_cell_[0] + 2), int(number_of_cells_[0])); ++i)
			for (size_t j = SMAX(int(body_upper_bound_cell_[1]) - 1, 0);
					j <= SMIN(int(body_upper_bound_cell_[1] + 1), int(number_of_cells_[1] - 1)); ++j) 
			{
				upper_bound_cells_.push_back(Vecu(i, j));
			}
	}
//=================================================================================================//
	void BoundingInYDirection::exec(Real dt)
	{
		
		//check lower bound
		for (size_t i = 0; i != lower_bound_cells_.size(); ++i) 
		{
			ListDataVector &list_data
				= cell_linked_lists_[lower_bound_cells_[i][0]][lower_bound_cells_[i][1]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				CheckLowerBound(list_data[num].first, list_data[num].second, dt);
		}

		//check upper bound
		for (size_t i = 0; i != upper_bound_cells_.size(); ++i) 
		{
			ListDataVector &list_data
				= cell_linked_lists_[upper_bound_cells_[i][0]][upper_bound_cells_[i][1]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				CheckUpperBound(list_data[num].first, list_data[num].second, dt);
		}
		
	}
//=================================================================================================//
	void BoundingInYDirection::parallel_exec(Real dt)
	{
		
		//check lower bound
		parallel_for(blocked_range<size_t>(0, lower_bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				ListDataVector &list_data
					= cell_linked_lists_[lower_bound_cells_[i][0]][lower_bound_cells_[i][1]]
					.particle_data_lists_;
				for (size_t num = 0; num < list_data.size(); ++num)
					CheckLowerBound(list_data[num].first, list_data[num].second, dt);
			}
		}, ap);

		//check upper bound
		parallel_for(blocked_range<size_t>(0, upper_bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				ListDataVector &list_data
					= cell_linked_lists_[upper_bound_cells_[i][0]][upper_bound_cells_[i][1]]
					.particle_data_lists_;
				for (size_t num = 0; num < list_data.size(); ++num)
					CheckUpperBound(list_data[num].first, list_data[num].second, dt);
			}
		}, ap);
		
	}
//=================================================================================================//
	BoundingInZDirection::BoundingInZDirection(SPHBody* body)
		: BoundingBodyDomain(body)
	{

	}
//=================================================================================================//
	void BoundingInZDirection::exec(Real dt)
	{

	}
//=================================================================================================//
	void BoundingInZDirection::parallel_exec(Real dt)
	{

	}
//=================================================================================================//
	void PeriodicConditionInXDirection::CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt)
	{
		Vecd particle_position = pnt;
		if (particle_position[0] > body_lower_bound_[0]
			&& particle_position[0] < (body_lower_bound_[0] + cell_spacing_))
		{
			Vecd translated_position = particle_position + periodic_translation_;
			Vecu cellpos = mesh_cell_linked_list_
				->GridIndexesFromPosition(translated_position);
			cell_linked_lists_[cellpos[0]][cellpos[1]].particle_data_lists_
				.push_back(make_pair(index_particle_i, translated_position));
		}
	}
//=================================================================================================//
	void PeriodicConditionInXDirection::CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt)
	{
		Vecd particle_position = pnt;
		if (particle_position[0] < body_upper_bound_[0]
			&& particle_position[0] > (body_upper_bound_[0] - cell_spacing_))
		{

			Vecd translated_position = particle_position - periodic_translation_;
			Vecu cellpos = mesh_cell_linked_list_
				->GridIndexesFromPosition(translated_position);
			cell_linked_lists_[cellpos[0]][cellpos[1]].particle_data_lists_
				.push_back(make_pair(index_particle_i, translated_position));
		}
	}
//=================================================================================================//
	void PeriodicConditionInYDirection::CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt)
	{
		Vecd particle_position = pnt;
		if (particle_position[1] > body_lower_bound_[1]
			&& particle_position[1] < (body_lower_bound_[1] + cell_spacing_))
		{
			Vecd translated_position = particle_position + periodic_translation_;
			Vecu cellpos = mesh_cell_linked_list_
				->GridIndexesFromPosition(translated_position);
			cell_linked_lists_[cellpos[0]][cellpos[1]].particle_data_lists_
				.push_back(make_pair(index_particle_i, translated_position));
		}
	}
//=================================================================================================//
	void PeriodicConditionInYDirection::CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt)
	{
		Vecd particle_position = pnt;
		if (particle_position[1] < body_upper_bound_[1]
			&& particle_position[1] > (body_upper_bound_[1] - cell_spacing_))
		{

			Vecd translated_position = particle_position - periodic_translation_;
			Vecu cellpos = mesh_cell_linked_list_
				->GridIndexesFromPosition(translated_position);
			cell_linked_lists_[cellpos[0]][cellpos[1]].particle_data_lists_
				.push_back(make_pair(index_particle_i, translated_position));
		}
	}
//=================================================================================================//
	void PeriodicConditionInZDirection
		::CheckLowerBound(size_t index_particle_i, Vecd pnt, Real dt)
	{

	}
//=================================================================================================//
	void PeriodicConditionInZDirection
		::CheckUpperBound(size_t index_particle_i, Vecd pnt, Real dt)
	{

	}
//=================================================================================================//
}