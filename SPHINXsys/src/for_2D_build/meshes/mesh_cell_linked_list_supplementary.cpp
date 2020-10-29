/**
 * @file 	mesh_cell_linked_list.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "mesh_cell_linked_list.h"
#include "base_kernel.h"
#include "base_body.h"
#include "base_particles.h"
#include "neighbor_relation.h"

namespace SPH {
	//=================================================================================================//
	CellList::CellList()
	{
		concurrent_particle_indexes_.reserve(12);
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList
		::ClearCellLists(Vecu& number_of_cells, matrix_cell cell_linked_lists)
	{
		parallel_for(blocked_range2d<size_t>(0, number_of_cells[0], 0, number_of_cells[1]),
			[&](const blocked_range2d<size_t>& r) {
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
						cell_linked_lists[i][j].concurrent_particle_indexes_.clear();
						cell_linked_lists[i][j].real_particle_indexes_.clear();
					}
			}, ap);
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList::UpdateSplitCellLists(SplitCellLists& split_cell_lists,
		Vecu& number_of_cells, matrix_cell cell_linked_lists)
	{
		//clear the data
		ClearSplitCellLists(split_cell_lists);

		parallel_for(blocked_range2d<size_t>(0, number_of_cells[0], 0, number_of_cells[1]),
			[&](const blocked_range2d<size_t>& r) {
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
						CellList& cell_list = cell_linked_lists[i][j];
						size_t real_particles_in_cell = cell_list.concurrent_particle_indexes_.size();
						if (real_particles_in_cell != 0) {
							for (size_t s = 0; s != real_particles_in_cell; ++s)
								cell_list.real_particle_indexes_.push_back(cell_list.concurrent_particle_indexes_[s]);
							split_cell_lists[transferMeshIndexTo1D(Vecu(3), Vecu(i % 3, j % 3))].push_back(&cell_linked_lists[i][j]);
						}
					}
			}, ap);
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList::UpdateCellListData(matrix_cell cell_linked_lists)
	{
		StdLargeVec<Vecd>& pos_n = base_particles_->pos_n_;
		parallel_for(blocked_range2d<size_t>(0, number_of_cells_[0], 0, number_of_cells_[1]),
			[&](const blocked_range2d<size_t>& r) {
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
						CellList& cell_list = cell_linked_lists[i][j];
						cell_list.cell_list_data_.clear();
						for (size_t s = 0; s != cell_list.concurrent_particle_indexes_.size(); ++s) {
							size_t particle_index = cell_list.concurrent_particle_indexes_[s];
							cell_list.cell_list_data_.emplace_back(make_pair(particle_index, pos_n[particle_index]));
						}
					}
			}, ap);
	}
	//=================================================================================================//
	CellList* MeshCellLinkedList::CellListFromIndex(Vecu cell_index)
	{
		return &cell_linked_lists_[cell_index[0]][cell_index[1]];
	}
	//=================================================================================================//
	void MeshCellLinkedList::allocateMeshDataMatrix()
	{
		Allocate2dArray(cell_linked_lists_, number_of_cells_);
	}
	//=================================================================================================//
	void MeshCellLinkedList::deleteMeshDataMatrix()
	{
		Delete2dArray(cell_linked_lists_, number_of_cells_);
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList
		::InsertACellLinkedListEntryAtALevel(size_t particle_index, Vecd& position, Vecu& cell_index, size_t level)
	{
		cell_linked_lists_levels_[level][cell_index[0]][cell_index[1]].concurrent_particle_indexes_
			.push_back(particle_index);
	}
	//=================================================================================================//
	CellList* MultilevelMeshCellLinkedList::CellListFromIndex(Vecu cell_index)
	{
		return &cell_linked_lists_levels_[0][cell_index[0]][cell_index[1]];
	}
	//=================================================================================================//
	void MeshCellLinkedList
		::InsertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position)
	{
		Vecu cellpos = GridIndexFromPosition(particle_position);
		cell_linked_lists_[cellpos[0]][cellpos[1]].concurrent_particle_indexes_.emplace_back(particle_index);
	}
	//=================================================================================================//
	void MeshCellLinkedList
		::InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position)
	{
		Vecu cellpos = GridIndexFromPosition(particle_position);
		cell_linked_lists_[cellpos[0]][cellpos[1]].cell_list_data_
			.emplace_back(make_pair(particle_index, particle_position));
	}
	//=================================================================================================//
	ListData MeshCellLinkedList::findNearestListDataEntry(Vecd& position)
	{
		Real min_distance = Infinity;
		ListData nearest_entry = std::make_pair(MaxSize_t, Vecd(Infinity));

		Vecu cell_location = GridIndexFromPosition(position);
		int i = (int)cell_location[0];
		int j = (int)cell_location[1];

		for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
		{
			for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
			{
				CellListDataVector& target_particles = cell_linked_lists_[l][m].cell_list_data_;
				for (size_t n = 0; n != target_particles.size(); ++n)
				{
					Real distance = (position - target_particles[n].second).norm();
					if (distance < min_distance)
					{
						min_distance = distance;
						nearest_entry = target_particles[n];
					}
				}
			}
		}
		return nearest_entry;
	}
	//=================================================================================================//
}
