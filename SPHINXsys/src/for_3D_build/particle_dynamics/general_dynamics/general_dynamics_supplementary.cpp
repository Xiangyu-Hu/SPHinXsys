/**
 * @file 	general_dynamics_supplementary.cpp
 * @brief 	supplementary functions declared in general_dynamics.h.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 */
#include "general_dynamics.h"

namespace SPH 
{
	//=================================================================================================//
	BoundingInAxisDirection::BoundingInAxisDirection(SPHBody* body, int axis_direction)
		: BoundingBodyDomain(body), axis_(axis_direction), second_axis_(SecondAxis(axis_direction)),
		third_axis_(ThirdAxis(axis_direction)) 	{}
	//=================================================================================================//
	PeriodicConditionInAxisDirection::PeriodicConditionInAxisDirection(SPHBody* body, int axis_direction) :
		BoundingInAxisDirection(body, axis_direction)
	{
		//allocate memory
		bound_cells_.resize(2);

		//lower bound cells
		for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_]) - 1, 0);
			k < (size_t)SMIN(int(body_upper_bound_cell_[third_axis_] + 2), int(number_of_cells_[third_axis_])); ++k)
		{

			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_]) - 1, 0);
				j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis_] + 2), int(number_of_cells_[second_axis_])); ++j)
			{

				for (size_t i = SMAX(int(body_lower_bound_cell_[axis_]) - 1, 0);
					i <= (size_t)SMIN(int(body_lower_bound_cell_[axis_] + 1), int(number_of_cells_[axis_] - 1)); ++i)
				{
					Vecu cell_position(0);
					cell_position[axis_] = i;
					cell_position[second_axis_] = j;
					cell_position[third_axis_] = k;
					bound_cells_[0].push_back(Vecu(cell_position));
				}
			}
		}

		//upper bound cells
		for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_]) - 1, 0);
			k < (size_t)SMIN(int(body_upper_bound_cell_[third_axis_] + 2), int(number_of_cells_[third_axis_])); ++k)
		{

			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_]) - 1, 0);
				j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis_] + 2), int(number_of_cells_[second_axis_])); ++j)
			{

				for (size_t i = SMAX(int(body_upper_bound_cell_[axis_]) - 1, 0);
					i <= (size_t)SMIN(int(body_upper_bound_cell_[axis_] + 1), int(number_of_cells_[axis_] - 1)); ++i)
				{
					Vecu cell_position(0);
					cell_position[axis_] = i;
					cell_position[second_axis_] = j;
					cell_position[third_axis_] = k;
					bound_cells_[1].push_back(Vecu(cell_position));
				}
			}
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicBounding::exec(Real dt)
	{
		//check lower bound
		CellVector& lower_bound_cells = bound_cells_[0];
		for (size_t i = 0; i != lower_bound_cells.size(); ++i) {
			IndexVector& particle_indexes
				= cell_linked_lists_[lower_bound_cells[i][0]][lower_bound_cells[i][1]][lower_bound_cells[i][2]]
				.real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num)
				checkLowerBound(particle_indexes[num], dt);
		}

		//check upper bound
		CellVector& upper_bound_cells = bound_cells_[1];
		for (size_t i = 0; i != upper_bound_cells.size(); ++i) {
			IndexVector& particle_indexes
				= cell_linked_lists_[upper_bound_cells[i][0]][upper_bound_cells[i][1]][upper_bound_cells[i][2]]
				.real_particle_indexes_;
			for (size_t num = 0; num < particle_indexes.size(); ++num)
				checkUpperBound(particle_indexes[num], dt);
		}
	}
	//=================================================================================================//
	void PeriodicConditionInAxisDirection::PeriodicBounding::parallel_exec(Real dt)
	{
		//check lower bound
		CellVector& lower_bound_cells = bound_cells_[0];
		parallel_for(blocked_range<size_t>(0, lower_bound_cells.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					IndexVector& particle_indexes
						= cell_linked_lists_[lower_bound_cells[i][0]][lower_bound_cells[i][1]][lower_bound_cells[i][2]]
						.real_particle_indexes_;
					for (size_t num = 0; num < particle_indexes.size(); ++num)
						checkLowerBound(particle_indexes[num], dt);
				}
			}, ap);

		//check upper bound
		CellVector& upper_bound_cells = bound_cells_[1];
		parallel_for(blocked_range<size_t>(0, upper_bound_cells.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					IndexVector& particle_indexes
						= cell_linked_lists_[upper_bound_cells[i][0]][upper_bound_cells[i][1]][upper_bound_cells[i][2]]
						.real_particle_indexes_;
					for (size_t num = 0; num < particle_indexes.size(); ++num)
						checkUpperBound(particle_indexes[num], dt);
				}
			}, ap);
	}
	//=================================================================================================//
	MirrorBoundaryConditionInAxisDirection
		::MirrorBoundaryConditionInAxisDirection(SPHBody* body, int axis_direction, bool positive)
		: BoundingInAxisDirection(body, axis_direction),
		bounding_(this->bound_cells_, body, axis_direction, positive),
		creating_ghost_particles_(this->ghost_particles_, this->bound_cells_, body, axis_direction, positive),
		updating_ghost_states_(this->ghost_particles_, this->bound_cells_, body, axis_direction, positive)
	{
		if (positive) {
			//upper bound cells
			for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_]) - 1, 0);
				k < (size_t)SMIN(int(body_upper_bound_cell_[third_axis_] + 2), int(number_of_cells_[third_axis_])); ++k)
			{

				for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_]) - 1, 0);
					j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis_] + 2), int(number_of_cells_[second_axis_])); ++j)
				{

					for (size_t i = SMAX(int(body_upper_bound_cell_[axis_]) - 1, 0);
						i <= (size_t)SMIN(int(body_upper_bound_cell_[axis_] + 1), int(number_of_cells_[axis_] - 1)); ++i)
					{
						Vecu cell_position(0);
						cell_position[axis_] = i;
						cell_position[second_axis_] = j;
						cell_position[third_axis_] = k;
						bound_cells_.push_back(Vecu(cell_position));
					}
				}
			}
		}
		else {
			//lower bound cells
			for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_]) - 1, 0);
				k < (size_t)SMIN(int(body_upper_bound_cell_[third_axis_] + 2), int(number_of_cells_[third_axis_])); ++k)
			{

				for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_]) - 1, 0);
					j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis_] + 2), int(number_of_cells_[second_axis_])); ++j)
				{

					for (size_t i = SMAX(int(body_lower_bound_cell_[axis_]) - 1, 0);
						i <= (size_t)SMIN(int(body_lower_bound_cell_[axis_] + 1), int(number_of_cells_[axis_] - 1)); ++i)
					{
						Vecu cell_position(0);
						cell_position[axis_] = i;
						cell_position[second_axis_] = j;
						cell_position[third_axis_] = k;
						bound_cells_.push_back(Vecu(cell_position));
					}
				}
			}
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::MirrorBounding
		::exec(Real dt)
	{
		for (size_t i = 0; i != bound_cells_.size(); ++i) {
			CellListDataVector& list_data
				= cell_linked_lists_[bound_cells_[i][0]][bound_cells_[i][1]][bound_cells_[i][2]]
				.cell_list_data_;
			for (size_t num = 0; num < list_data.size(); ++num)
				checking_bound_(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void MirrorBoundaryConditionInAxisDirection::MirrorBounding
		::parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					CellListDataVector& list_data
						= cell_linked_lists_[bound_cells_[i][0]][bound_cells_[i][1]][bound_cells_[i][2]]
						.cell_list_data_;
					for (size_t num = 0; num < list_data.size(); ++num)
						checking_bound_(list_data[num].first, dt);
				}
			}, ap);
	}
	//=================================================================================================//
}