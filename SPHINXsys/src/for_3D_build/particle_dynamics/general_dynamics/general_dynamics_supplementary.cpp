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
		third_axis_(ThirdAxis(axis_direction))
	{
		//lower bound cells
		for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_] - 1), 0);
			k < SMIN(int(body_upper_bound_cell_[third_axis_] + 2), 
				int(number_of_cells_[third_axis_])); ++k)
		{

			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_] - 1), 0);
				j < SMIN(int(body_upper_bound_cell_[second_axis_] + 2), 
					int(number_of_cells_[second_axis_])); ++j)
			{

				for (size_t i = SMAX(int(body_lower_bound_cell_[axis_]) - 1, 0);
					i <= SMIN(int(body_lower_bound_cell_[axis_] + 1), int(number_of_cells_[axis_] - 1)); ++i)
				{
					Vecu cell_position(0);
					cell_position[axis_] = i;
					cell_position[second_axis_] = j;
					cell_position[third_axis_] = k;
					lower_bound_cells_.push_back(Vecu(cell_position));
				}
			}
		}

		//upper bound cells
		for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_] - 1), 0);
			k < SMIN(int(body_upper_bound_cell_[third_axis_] + 2), 
				int(number_of_cells_[third_axis_])); ++k)
		{

			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_] - 1), 0);
				j < SMIN(int(body_upper_bound_cell_[second_axis_] + 2), 
					int(number_of_cells_[second_axis_])); ++j)
			{

				for (size_t i = SMAX(int(body_upper_bound_cell_[axis_]) - 1, 0);
					i <= SMIN(int(body_upper_bound_cell_[axis_] + 1), 
						int(number_of_cells_[axis_] - 1)); ++i)
				{
					Vecu cell_position(0);
					cell_position[axis_] = i;
					cell_position[second_axis_] = j;
					cell_position[third_axis_] = k;
					upper_bound_cells_.push_back(Vecu(cell_position));
				}
			}
		}
	}
	//=================================================================================================//
	void BoundingInAxisDirection::exec(Real dt)
	{
		//check lower bound
		for (size_t i = 0; i != lower_bound_cells_.size(); ++i) {
			ConcurrentListDataVector& list_data
				= cell_linked_lists_[lower_bound_cells_[i][0]][lower_bound_cells_[i][1]][lower_bound_cells_[i][2]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				CheckLowerBound(list_data[num].first, list_data[num].second, dt);
		}

		//check upper bound
		for (size_t i = 0; i != upper_bound_cells_.size(); ++i) {
			ConcurrentListDataVector& list_data
				= cell_linked_lists_[upper_bound_cells_[i][0]][upper_bound_cells_[i][1]][upper_bound_cells_[i][2]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				CheckUpperBound(list_data[num].first, list_data[num].second, dt);
		}
	}
	//=================================================================================================//
	void BoundingInAxisDirection::parallel_exec(Real dt)
	{
		//check lower bound
		parallel_for(blocked_range<size_t>(0, lower_bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					ConcurrentListDataVector& list_data
						= cell_linked_lists_[lower_bound_cells_[i][0]][lower_bound_cells_[i][1]][lower_bound_cells_[i][2]]
						.particle_data_lists_;
					for (size_t num = 0; num < list_data.size(); ++num)
						CheckLowerBound(list_data[num].first, list_data[num].second, dt);
				}
			}, ap);

		//check upper bound
		parallel_for(blocked_range<size_t>(0, upper_bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					ConcurrentListDataVector& list_data
						= cell_linked_lists_[upper_bound_cells_[i][0]][upper_bound_cells_[i][1]][upper_bound_cells_[i][2]]
						.particle_data_lists_;
					for (size_t num = 0; num < list_data.size(); ++num)
						CheckUpperBound(list_data[num].first, list_data[num].second, dt);
				}
			}, ap);
	}
	//=================================================================================================//
	MirrorBoundingInAxisDirection
		::MirrorBoundingInAxisDirection(SPHBody* body, int axis_direction, bool positive)
		: BoundingBodyDomain(body), axis_(axis_direction), second_axis_(SecondAxis(axis_direction)),
		third_axis_(ThirdAxis(axis_direction))
	{
		if (positive) {
			//upper bound cells
			for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_] - 1), 0);
				k < SMIN(int(body_upper_bound_cell_[third_axis_] + 2),
					int(number_of_cells_[third_axis_])); ++k)
			{

				for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_] - 1), 0);
					j < SMIN(int(body_upper_bound_cell_[second_axis_] + 2),
						int(number_of_cells_[second_axis_])); ++j)
				{

					for (size_t i = SMAX(int(body_upper_bound_cell_[axis_]) - 1, 0);
						i <= SMIN(int(body_upper_bound_cell_[axis_] + 1),
							int(number_of_cells_[axis_] - 1)); ++i)
					{
						Vecu cell_position(0);
						cell_position[axis_] = i;
						cell_position[second_axis_] = j;
						cell_position[third_axis_] = k;
						bound_cells_.push_back(Vecu(cell_position));
					}
				}
			}
			checking_bound_ = std::bind(&MirrorBoundingInAxisDirection::CheckUpperBound, this, _1, _2);
		}
		else {
			//lower bound cells
			for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis_] - 1), 0);
				k < SMIN(int(body_upper_bound_cell_[third_axis_] + 2),
					int(number_of_cells_[third_axis_])); ++k)
			{

				for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis_] - 1), 0);
					j < SMIN(int(body_upper_bound_cell_[second_axis_] + 2),
						int(number_of_cells_[second_axis_])); ++j)
				{

					for (size_t i = SMAX(int(body_lower_bound_cell_[axis_]) - 1, 0);
						i <= SMIN(int(body_lower_bound_cell_[axis_] + 1), int(number_of_cells_[axis_] - 1)); ++i)
					{
						Vecu cell_position(0);
						cell_position[axis_] = i;
						cell_position[second_axis_] = j;
						cell_position[third_axis_] = k;
						bound_cells_.push_back(Vecu(cell_position));
					}
				}
			}
			checking_bound_ = std::bind(&MirrorBoundingInAxisDirection::CheckLowerBound, this, _1, _2);
		}
	}
	//=================================================================================================//
	void MirrorBoundingInAxisDirection::exec(Real dt)
	{
		for (size_t i = 0; i != bound_cells_.size(); ++i) {
			ConcurrentListDataVector& list_data
				= cell_linked_lists_[bound_cells_[i][0]][bound_cells_[i][1]][bound_cells_[i][2]]
				.particle_data_lists_;
			for (size_t num = 0; num < list_data.size(); ++num)
				checking_bound_(list_data[num].first, dt);
		}
	}
	//=================================================================================================//
	void MirrorBoundingInAxisDirection::parallel_exec(Real dt)
	{
		parallel_for(blocked_range<size_t>(0, bound_cells_.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i < r.end(); ++i) {
					ConcurrentListDataVector& list_data
						= cell_linked_lists_[bound_cells_[i][0]][bound_cells_[i][1]][bound_cells_[i][2]]
						.particle_data_lists_;
					for (size_t num = 0; num < list_data.size(); ++num)
						checking_bound_(list_data[num].first, dt);
				}
			}, ap);
	}
	//=================================================================================================//
}