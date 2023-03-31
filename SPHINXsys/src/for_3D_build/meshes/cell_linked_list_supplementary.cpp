/**
 * @file 	cell_linked_list_supplementary.cpp
 * @author	Luhui Han, Chi Zhang and Xiangyu Hu
 */

#include "cell_linked_list.h"

#include "base_particles.hpp"
#include "mesh_iterators.hpp"

namespace SPH
{
	//=================================================================================================//
	void CellLinkedList ::allocateMeshDataMatrix()
	{
		Allocate3dArray(cell_index_lists_, all_cells_);
		Allocate3dArray(cell_data_lists_, all_cells_);

		mesh_parallel_for(MeshRange(Array3i::Zero(), all_cells_),
						  [&](int i, int j, int k)
						  {
							  cell_index_lists_[i][j][k].reserve(36);
							  cell_data_lists_[i][j][k].reserve(36);
						  });
	}
	//=================================================================================================//
	void CellLinkedList ::deleteMeshDataMatrix()
	{
		Delete3dArray(cell_index_lists_, all_cells_);
		Delete3dArray(cell_data_lists_, all_cells_);
	}
	//=================================================================================================//
	void CellLinkedList::clearCellLists()
	{
		mesh_parallel_for(MeshRange(Array3i::Zero(), all_cells_),
						  [&](int i, int j, int k)
						  {
							  cell_index_lists_[i][j][k].clear();
						  });
	}
	//=================================================================================================//
	void CellLinkedList::UpdateCellListData(BaseParticles &base_particles)
	{
		StdLargeVec<Vecd> &pos = base_particles.pos_;
		StdLargeVec<Real> &Vol = base_particles.Vol_;
		mesh_parallel_for(MeshRange(Array3i::Zero(), all_cells_),
						  [&](int i, int j, int k)
						  {
							  cell_data_lists_[i][j][k].clear();
							  ConcurrentIndexVector &cell_list = cell_index_lists_[i][j][k];
							  for (size_t s = 0; s != cell_list.size(); ++s)
							  {
								  size_t index = cell_list[s];
								  cell_data_lists_[i][j][k].emplace_back(std::make_tuple(index, pos[index], Vol[index]));
							  }
						  });
	}
	//=================================================================================================//
	void CellLinkedList::updateSplitCellLists(SplitCellLists &split_cell_lists)
	{
		clearSplitCellLists(split_cell_lists);
		mesh_parallel_for(MeshRange(Array3i::Zero(), all_cells_),
						  [&](int i, int j, int k)
						  {
							  size_t real_particles_in_cell = cell_index_lists_[i][j][k].size();
							  if (real_particles_in_cell != 0)
							  {
								  split_cell_lists[transferMeshIndexTo1D(Array3i(3, 3, 3), Array3i(i % 3, j % 3, k % 3))]
									  .push_back(&cell_index_lists_[i][j][k]);
							  }
						  });
	}
	//=================================================================================================//
	void CellLinkedList ::insertParticleIndex(size_t particle_index, const Vecd &particle_position)
	{
		Array3i cell_pos = CellIndexFromPosition(particle_position);
		cell_index_lists_[cell_pos[0]][cell_pos[1]][cell_pos[2]].emplace_back(particle_index);
	}
	//=================================================================================================//
	void CellLinkedList ::InsertListDataEntry(size_t particle_index, const Vecd &particle_position, Real volumetric)
	{
		Array3i cell_pos = CellIndexFromPosition(particle_position);
		cell_data_lists_[cell_pos[0]][cell_pos[1]][cell_pos[2]].emplace_back(
			std::make_tuple(particle_index, particle_position, volumetric));
	}
	//=================================================================================================//
	ListData CellLinkedList::findNearestListDataEntry(const Vecd &position)
	{
		Real min_distance = Infinity;
		ListData nearest_entry = std::make_tuple(MaxSize_t, Infinity * Vecd::Ones(), Infinity);

		Array3i cell_location = CellIndexFromPosition(position);
		int i = (int)cell_location[0];
		int j = (int)cell_location[1];
		int k = (int)cell_location[2];

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(all_cells_[0]) - 1); ++l)
		{
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(all_cells_[1]) - 1); ++m)
			{
				for (int q = SMAX(k - 1, 0); q <= SMIN(k + 1, int(all_cells_[2]) - 1); ++q)
				{
					ListDataVector &target_particles = cell_data_lists_[l][m][q];
					for (const ListData &list_data : target_particles)
					{
						Real distance = (position - std::get<1>(list_data)).norm();
						if (distance < min_distance)
						{
							min_distance = distance;
							nearest_entry = list_data;
						}
					}
				}
			}
		}
		return nearest_entry;
	}
	//=================================================================================================//
	void CellLinkedList::
		tagBodyPartByCell(ConcurrentCellLists &cell_lists, std::function<bool(Vecd, Real)> &check_included)
	{
		mesh_parallel_for(MeshRange(Array3i::Zero(), all_cells_),
						  [&](int i, int j, int k)
						  {
							  bool is_included = false;
							  mesh_for_each(
								  Array3i::Zero().max(Array3i(i, j, k) - Array3i::Ones()),
								  all_cells_.min(Array3i(i, j, k) + Array3i::Ones()),
								  [&](int l, int m, int n)
								  {
									  if (check_included(CellPositionFromIndex(Array3i(l, m, n)), grid_spacing_))
									  {
										  is_included = true;
									  }
								  });
							  if (is_included == true)
								  cell_lists.push_back(&cell_index_lists_[i][j][k]);
						  });
	}
	//=================================================================================================//
	void CellLinkedList::
		tagBoundingCells(StdVec<CellLists> &cell_data_lists, BoundingBox &bounding_bounds, int axis)
	{
		int second_axis = SecondAxis(axis);
		int third_axis = ThirdAxis(axis);
		Array3i body_lower_bound_cell_ = CellIndexFromPosition(bounding_bounds.first_);
		Array3i body_upper_bound_cell_ = CellIndexFromPosition(bounding_bounds.second_);

		// lower bound cells
		for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis]) - 1, 0);
			 k < (size_t)SMIN(int(body_upper_bound_cell_[third_axis] + 2), int(all_cells_[third_axis])); ++k)
		{

			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis]) - 1, 0);
				 j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis] + 2), int(all_cells_[second_axis])); ++j)
			{

				for (size_t i = SMAX(int(body_lower_bound_cell_[axis]) - 1, 0);
					 i <= (size_t)SMIN(int(body_lower_bound_cell_[axis] + 1), int(all_cells_[axis] - 1)); ++i)
				{
					Array3i cell_position = Array3i::Zero();
					cell_position[axis] = i;
					cell_position[second_axis] = j;
					cell_position[third_axis] = k;
					cell_data_lists[0].first.push_back(&cell_index_lists_[cell_position[0]][cell_position[1]][cell_position[2]]);
					cell_data_lists[0].second.push_back(&cell_data_lists_[cell_position[0]][cell_position[1]][cell_position[2]]);
				}
			}
		}

		// upper bound cells
		for (size_t k = SMAX(int(body_lower_bound_cell_[third_axis]) - 1, 0);
			 k < (size_t)SMIN(int(body_upper_bound_cell_[third_axis] + 2), int(all_cells_[third_axis])); ++k)
		{

			for (size_t j = SMAX(int(body_lower_bound_cell_[second_axis]) - 1, 0);
				 j < (size_t)SMIN(int(body_upper_bound_cell_[second_axis] + 2), int(all_cells_[second_axis])); ++j)
			{

				for (size_t i = SMAX(int(body_upper_bound_cell_[axis]) - 1, 0);
					 i <= (size_t)SMIN(int(body_upper_bound_cell_[axis] + 1), int(all_cells_[axis] - 1)); ++i)
				{
					Array3i cell_position = Array3i::Zero();
					cell_position[axis] = i;
					cell_position[second_axis] = j;
					cell_position[third_axis] = k;
					cell_data_lists[1].first.push_back(&cell_index_lists_[cell_position[0]][cell_position[1]][cell_position[2]]);
					cell_data_lists[1].second.push_back(&cell_data_lists_[cell_position[0]][cell_position[1]][cell_position[2]]);
				}
			}
		}
	}
	//=================================================================================================//
	void CellLinkedList::writeMeshFieldToPlt(std::ofstream &output_file)
	{
		Array3i number_of_operation = all_cells_;

		output_file << "\n";
		output_file << "title='View'"
					<< "\n";
		output_file << "variables= "
					<< "x, "
					<< "y, "
					<< "z, "
					<< "particles_in_cell "
					<< "\n";
		output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
					<< "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

		for (int k = 0; k != number_of_operation[2]; ++k)
			for (int j = 0; j != number_of_operation[1]; ++j)
			{
				for (int i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = CellPositionFromIndex(Array3i(i, j, k));
					output_file << data_position[0] << " ";
				}
				output_file << " \n";
			}

		for (int k = 0; k != number_of_operation[2]; ++k)
			for (int j = 0; j != number_of_operation[1]; ++j)
			{
				for (int i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = CellPositionFromIndex(Array3i(i, j, k));
					output_file << data_position[1] << " ";
				}
				output_file << " \n";
			}

		for (int k = 0; k != number_of_operation[2]; ++k)
			for (int j = 0; j != number_of_operation[1]; ++j)
			{
				for (int i = 0; i != number_of_operation[0]; ++i)
				{
					Vecd data_position = CellPositionFromIndex(Array3i(i, j, k));
					output_file << data_position[2] << " ";
				}
				output_file << " \n";
			}

		for (int k = 0; k != number_of_operation[2]; ++k)
			for (int j = 0; j != number_of_operation[1]; ++j)
			{
				for (int i = 0; i != number_of_operation[0]; ++i)
				{
					output_file << cell_index_lists_[i][j][k].size() << " ";
				}
				output_file << " \n";
			}
	}
	//=================================================================================================//
}
