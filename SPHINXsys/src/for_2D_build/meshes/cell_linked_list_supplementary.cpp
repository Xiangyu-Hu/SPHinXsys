/**
 * @file 	cell_linked_list.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "cell_linked_list.h"
#include "base_kernel.h"
#include "base_body.h"
#include "base_particles.h"
#include "neighbor_relation.h"

namespace SPH
{
	//=================================================================================================//
	CellList::CellList()
	{
		concurrent_particle_indexes_.reserve(12);
	}
	//=================================================================================================//
	void CellLinkedList::allocateMeshDataMatrix()
	{
		Allocate2dArray(cell_linked_lists_, number_of_cells_);
	}
	//=================================================================================================//
	void CellLinkedList::deleteMeshDataMatrix()
	{
		Delete2dArray(cell_linked_lists_, number_of_cells_);
	}
	//=================================================================================================//
	void CellLinkedList::clearCellLists()
	{
		parallel_for(
			blocked_range2d<size_t>(0, number_of_cells_[0], 0, number_of_cells_[1]),
			[&](const blocked_range2d<size_t> &r)
			{
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
					{
						cell_linked_lists_[i][j].concurrent_particle_indexes_.clear();
						cell_linked_lists_[i][j].real_particle_indexes_.clear();
					}
			},
			ap);
	}
	//=================================================================================================//
	void CellLinkedList::UpdateCellListData()
	{
		StdLargeVec<Vecd> &pos_n = base_particles_->pos_n_;
		parallel_for(
			blocked_range2d<size_t>(0, number_of_cells_[0], 0, number_of_cells_[1]),
			[&](const blocked_range2d<size_t> &r)
			{
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
					{
						CellList &cell_list = cell_linked_lists_[i][j];
						cell_list.cell_list_data_.clear();
						for (size_t s = 0; s != cell_list.concurrent_particle_indexes_.size(); ++s)
						{
							size_t particle_index = cell_list.concurrent_particle_indexes_[s];
							cell_list.cell_list_data_.emplace_back(std::make_pair(particle_index, pos_n[particle_index]));
						}
					}
			},
			ap);
	}
	//=================================================================================================//
	void CellLinkedList::updateSplitCellLists(SplitCellLists &split_cell_lists)
	{
		//clear the data
		clearSplitCellLists(split_cell_lists);

		parallel_for(
			blocked_range2d<size_t>(0, number_of_cells_[0], 0, number_of_cells_[1]),
			[&](const blocked_range2d<size_t> &r)
			{
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
					{
						CellList &cell_list = cell_linked_lists_[i][j];
						size_t real_particles_in_cell = cell_list.concurrent_particle_indexes_.size();
						if (real_particles_in_cell != 0)
						{
							for (size_t s = 0; s != real_particles_in_cell; ++s)
								cell_list.real_particle_indexes_.push_back(cell_list.concurrent_particle_indexes_[s]);
							split_cell_lists[transferMeshIndexTo1D(Vecu(3), Vecu(i % 3, j % 3))].push_back(&cell_linked_lists_[i][j]);
						}
					}
			},
			ap);
	}
	//=================================================================================================//
	void CellLinkedList ::insertACellLinkedParticleIndex(size_t particle_index, const Vecd &particle_position)
	{
		Vecu cellpos = CellIndexFromPosition(particle_position);
		cell_linked_lists_[cellpos[0]][cellpos[1]].concurrent_particle_indexes_.emplace_back(particle_index);
	}
	//=================================================================================================//
	void CellLinkedList ::InsertACellLinkedListDataEntry(size_t particle_index, const Vecd &particle_position)
	{
		Vecu cellpos = CellIndexFromPosition(particle_position);
		cell_linked_lists_[cellpos[0]][cellpos[1]].cell_list_data_.emplace_back(std::make_pair(particle_index, particle_position));
	}
	//=================================================================================================//
	ListData CellLinkedList::findNearestListDataEntry(const Vecd &position)
	{
		Real min_distance = Infinity;
		ListData nearest_entry = std::make_pair(MaxSize_t, Vecd(Infinity));

		Vecu cell_location = CellIndexFromPosition(position);
		int i = (int)cell_location[0];
		int j = (int)cell_location[1];

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
		{
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
			{
				ListDataVector &target_particles = cell_linked_lists_[l][m].cell_list_data_;
				for (const ListData &list_data : target_particles)
				{
					Real distance = (position - list_data.second).norm();
					if (distance < min_distance)
					{
						min_distance = distance;
						nearest_entry = list_data;
					}
				}
			}
		}
		return nearest_entry;
	}
	//=================================================================================================//
	void CellLinkedList::
		tagBodyPartByCell(CellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included)
	{
		for (int i = 0; i < (int)number_of_cells_[0]; ++i)
			for (int j = 0; j < (int)number_of_cells_[1]; ++j)
			{
				bool is_included = false;
				for (int k = SMAX(i - 1, 0); k <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++k)
					for (int l = SMAX(j - 1, 0); l <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++l)
					{
						if (check_included(CellPositionFromIndex(Vecu(k, l)), grid_spacing_))
						{
							is_included = true;
						}
					}
				if (is_included == true)
					cell_lists.push_back(&cell_linked_lists_[i][j]);
			}
	}
	//=================================================================================================//
	void CellLinkedList::
		tagBodyDomainBoundingCells(StdVec<CellLists> &cell_lists, BoundingBox &body_domain_bounds, int axis)
	{
		int second_axis = SecondAxis(axis);
		Vecu body_lower_bound_cell_ = CellIndexFromPosition(body_domain_bounds.first);
		Vecu body_upper_bound_cell_ = CellIndexFromPosition(body_domain_bounds.second);

		//lower bound cells
		for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis]) - 1, 0);
			 j <= (size_t)SMIN(int(body_upper_bound_cell_[second_axis] + 1), int(number_of_cells_[second_axis] - 1)); ++j)
			for (size_t i = SMAX(int(body_lower_bound_cell_[axis]) - 1, 0);
				 i <= (size_t)SMIN(int(body_lower_bound_cell_[axis] + 1), int(number_of_cells_[axis] - 1)); ++i)
			{
				Vecu cell_position(0);
				cell_position[axis] = i;
				cell_position[second_axis] = j;
				cell_lists[0].push_back(&cell_linked_lists_[cell_position[0]][cell_position[1]]);
			}

		//upper bound cells
		for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis]) - 1, 0);
			 j <= (size_t)SMIN(int(body_upper_bound_cell_[second_axis] + 1), int(number_of_cells_[second_axis] - 1)); ++j)
			for (size_t i = SMAX(int(body_upper_bound_cell_[axis]) - 1, 0);
				 i <= (size_t)SMIN(int(body_upper_bound_cell_[axis] + 1), int(number_of_cells_[axis] - 1)); ++i)
			{
				Vecu cell_position(0);
				cell_position[axis] = i;
				cell_position[second_axis] = j;
				cell_lists[1].push_back(&cell_linked_lists_[cell_position[0]][cell_position[1]]);
			}
	}
	//=================================================================================================//
	void CellLinkedList::
		tagMirrorBoundingCells(CellLists &cell_lists, BoundingBox &body_domain_bounds, int axis, bool positive)
	{
		int second_axis = SecondAxis(axis);
		Vecu body_lower_bound_cell_ = CellIndexFromPosition(body_domain_bounds.first);
		Vecu body_upper_bound_cell_ = CellIndexFromPosition(body_domain_bounds.second);

		if (positive)
		{
			//upper bound cells
			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis]) - 1, 0);
				 j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis] + 2), int(number_of_cells_[second_axis])); ++j)
				for (size_t i = SMAX(int(body_upper_bound_cell_[axis]) - 1, 0);
					 i <= (size_t)SMIN(int(body_upper_bound_cell_[axis] + 1), int(number_of_cells_[axis] - 1)); ++i)
				{
					Vecu cell_position(0);
					cell_position[axis] = i;
					cell_position[second_axis] = j;
					cell_lists.push_back(&cell_linked_lists_[cell_position[0]][cell_position[1]]);
				}
		}
		else
		{
			//lower bound cells
			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis]) - 1, 0);
				 j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis] + 2), int(number_of_cells_[second_axis])); ++j)
				for (size_t i = SMAX(int(body_lower_bound_cell_[axis]) - 1, 0);
					 i <= (size_t)SMIN(int(body_lower_bound_cell_[axis] + 1), int(number_of_cells_[axis] - 1)); ++i)
				{
					Vecu cell_position(0);
					cell_position[axis] = i;
					cell_position[second_axis] = j;
					cell_lists.push_back(&cell_linked_lists_[cell_position[0]][cell_position[1]]);
				}
		}
	}
	//=============================================================================================//
	void CellLinkedList::writeMeshFieldToPlt(std::ofstream &output_file)
	{
		Vecu number_of_operation = number_of_cells_;

		output_file << "\n";
		output_file << "title='View'"
					<< "\n";
		output_file << "variables= "
					<< "x, "
					<< "y, "
					<< "particles_in_cell "
					<< "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << 1
					<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd data_position = CellPositionFromIndex(Vecu(i, j));
				output_file << data_position[0] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				Vecd data_position = CellPositionFromIndex(Vecu(i, j));
				output_file << data_position[1] << " ";
			}
			output_file << " \n";
		}

		for (size_t j = 0; j != number_of_operation[1]; ++j)
		{
			for (size_t i = 0; i != number_of_operation[0]; ++i)
			{
				output_file << cell_linked_lists_[i][j].concurrent_particle_indexes_.size() << " ";
			}
			output_file << " \n";
		}
	}
	//=================================================================================================//
}
