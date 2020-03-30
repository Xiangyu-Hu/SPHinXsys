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
	CellList::CellList() : cell_location_(0)
	{
		particle_data_lists_.reserve(12);
	}
	//=================================================================================================//
	void BaseMeshCellLinkedList
		::ClearCellLists(Vecu& number_of_cells, matrix_cell cell_linked_lists)
	{
		parallel_for(blocked_range2d<size_t>(0, number_of_cells[0], 0, number_of_cells[1]),
			[&](const blocked_range2d<size_t>& r) {
				for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
					for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
						cell_linked_lists[i][j].particle_data_lists_.clear();
						cell_linked_lists[i][j].real_particle_count_ = 0;
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

		parallel_for(blocked_range<size_t>(0, split_cell_lists.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num) {
					/** convert 1d vector index to mesh index. */
					Vec2u split_mesh_index = transfer1DtoMeshIndex(Vec2u(3, 3), num);
					size_t l = split_mesh_index[0];
					size_t m = split_mesh_index[1];

					Vecu starting_index(2 - l, 2 - m);
					Vecu number_of_operation = (number_of_cells + starting_index) / 3;
					for (size_t i = 0; i != number_of_operation[0]; ++i)
						for (size_t j = 0; j != number_of_operation[1]; ++j) {
							CellList& cell_list = cell_linked_lists[3 * i + l][3 * j + m];
							size_t real_particles_in_cell = cell_list.particle_data_lists_.size();
							if (real_particles_in_cell != 0) {
								cell_list.real_particle_count_ = real_particles_in_cell;
								for (int s = 0; s != real_particles_in_cell; ++s)
									cell_list.real_particle_indexes_.push_back(cell_list.particle_data_lists_[s].first);
								split_cell_lists[num].push_back(&cell_linked_lists[3 * i + l][3 * j + m]);
							}
						}
				}
			}, ap);
	}
	//=================================================================================================//
	CellList* MeshCellLinkedList::getCellList(Vecu cell_index)
	{
		return &cell_linked_lists_[cell_index[0]][cell_index[1]];
	}
	//=================================================================================================//
	void MeshCellLinkedList::AllocateMeshDataMatrix()
	{
		Allocate2dArray(cell_linked_lists_, number_of_cells_);
		for (size_t i = 0; i != number_of_cells_[0]; ++i)
			for (size_t j = 0; j != number_of_cells_[1]; ++j) {
				cell_linked_lists_[i][j].setCellInformation(Vecu(i,j));
			}
	}
	//=================================================================================================//
	void MeshCellLinkedList::DeleteMeshDataMatrix()
	{
		Delete2dArray(cell_linked_lists_, number_of_cells_);
	}
	//=================================================================================================//
	void MeshCellLinkedList::UpdateInnerConfiguration(ParticleConfiguration& inner_configuration)
	{
		StdLargeVec<BaseParticleData> &base_particle_data 
			= body_->base_particles_->base_particle_data_;

		parallel_for(blocked_range<size_t>(0, body_->number_of_particles_),
			[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num) {
					Vecu cell_location 
						= GridIndexesFromPosition(base_particle_data[num].pos_n_);
					int i = (int)cell_location[0];
					int j = (int)cell_location[1];

					Neighborhood& neighborhood = inner_configuration[num];
					NeighborList& neighbor_list = std::get<0>(neighborhood);
					size_t previous_count_of_neigbors = std::get<2>(neighborhood);

					for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
						for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
						{
							ConcurrentListDataVector& target_particles = cell_linked_lists_[l][m].particle_data_lists_;
							for (size_t n = 0; n != target_particles.size(); ++n)
							{
								//displacement pointing from neighboring particle to origin particle
								Vecd displacement = base_particle_data[num].pos_n_ - target_particles[n].second;
								if (displacement.norm() <= cutoff_radius_ && num != target_particles[n].first)
								{
									std::get<1>(neighborhood) >= previous_count_of_neigbors ?
										neighbor_list.push_back(new NeighborRelation(base_particle_data, *kernel_,
											displacement, num, target_particles[n].first))
										: neighbor_list[std::get<1>(neighborhood)]->resetRelation(base_particle_data,
											*kernel_, displacement, num, target_particles[n].first);
									std::get<1>(neighborhood)++;
								}
							}
						}
					size_t current_count_of_neighbors = std::get<1>(neighborhood);
					neighbor_list.resize(current_count_of_neighbors);
					std::get<2>(neighborhood) = current_count_of_neighbors;
					std::get<1>(neighborhood) = 0;
			}
		}, ap);
	}
	//=================================================================================================//
	void MeshCellLinkedList::UpdateInteractionConfiguration(SPHBodyVector interacting_bodies)
	{
		StdLargeVec<BaseParticleData> &base_particle_data = body_->base_particles_->base_particle_data_;
		ContatcParticleConfiguration& current_contact_configuration = body_->contact_configuration_;
		ContactParticles& indexes_contact_particles = body_->indexes_contact_particles_;
		SPHBodyVector contact_bodies = body_->contact_map_.second;
		IndexVector contact_configuration_index(interacting_bodies.size());

		for (size_t intertaction_body_num = 0;
			intertaction_body_num < interacting_bodies.size(); ++intertaction_body_num) {
			//clear previous interacting particles
			//build up contact configuration indexes
			for (size_t contact_body_num = 0;
				contact_body_num < contact_bodies.size(); ++contact_body_num) {
				if (interacting_bodies[intertaction_body_num] == contact_bodies[contact_body_num]) {
					contact_configuration_index[intertaction_body_num] = contact_body_num;
					indexes_contact_particles[contact_body_num].clear();
				}
			}

			BaseMeshCellLinkedList &target_mesh_cell_linked_list
				= *(interacting_bodies[intertaction_body_num]->base_mesh_cell_linked_list_);
			Vecu target_number_of_cells = target_mesh_cell_linked_list.getNumberOfCells();
			int search_range
				= ComputingSearchRage(body_->refinement_level_,
					interacting_bodies[intertaction_body_num]->refinement_level_);
			Kernel &current_kernel = ChoosingKernel(body_->kernel_,
				interacting_bodies[intertaction_body_num]->kernel_);
			Real cutoff_radius = current_kernel.GetCutOffRadius();

			parallel_for(blocked_range<size_t>(0, body_->number_of_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num) {

					Vecu target_cell_index = target_mesh_cell_linked_list
						.GridIndexesFromPosition(base_particle_data[num].pos_n_);
					int i = (int)target_cell_index[0];
					int j = (int)target_cell_index[1];

					matrix_cell target_cell_linked_lists
						= target_mesh_cell_linked_list.getCellLinkedLists();
					size_t contact_body_num
						= contact_configuration_index[intertaction_body_num];

					Neighborhood& neighborhood = current_contact_configuration[contact_body_num][num];
					NeighborList& neighbor_list = std::get<0>(neighborhood);
					size_t previous_count_of_neigbors = std::get<2>(neighborhood);

					for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
						{
							ConcurrentListDataVector &target_particles
								= target_cell_linked_lists[l][m].particle_data_lists_;
							for (size_t n = 0; n < target_particles.size(); n++)
							{
								//displacement pointing from neighboring particle to origin particle
								Vecd displacement = base_particle_data[num].pos_n_
									- target_particles[n].second;
								if (displacement.norm() <= cutoff_radius)
								{
									std::get<1>(neighborhood) >= previous_count_of_neigbors ?
										neighbor_list.push_back(new NeighborRelation(base_particle_data, current_kernel,
											displacement, num, target_particles[n].first))
										: neighbor_list[std::get<1>(neighborhood)]->resetRelation(base_particle_data,
											current_kernel, displacement, num, target_particles[n].first);
									std::get<1>(neighborhood)++;
								}
							}
						}
					size_t current_count_of_neighbors = std::get<1>(neighborhood);
					neighbor_list.resize(current_count_of_neighbors);
					std::get<2>(neighborhood) = current_count_of_neighbors;
					std::get<1>(neighborhood) = 0;
					if (current_count_of_neighbors != 0)
						indexes_contact_particles[contact_body_num].push_back(num);
				}
			}, ap);
		}
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList
		::InsertACellLinkedListEntryAtALevel(size_t particle_index, Vecd& position, Vecu& cell_index, size_t level)
	{
		cell_linked_lists_levels_[level][cell_index[0]][cell_index[1]].particle_data_lists_
			.push_back(make_pair(particle_index, position));
	}
	//=================================================================================================//
	CellList* MultilevelMeshCellLinkedList::getCellList(Vecu cell_index)
	{
		return &cell_linked_lists_levels_[0][cell_index[0]][cell_index[1]];
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList
		::UpdateInnerConfiguration(ParticleConfiguration& inner_configuration)
	{
		StdLargeVec<BaseParticleData>& base_particle_data = base_particles_->base_particle_data_;
		for(size_t level = 0; level != total_levels_; ++level) {
			matrix_cell cell_linked_lists = cell_linked_lists_levels_[level];
			Real cell_spacing = cell_spacing_levels_[level];

			SplitCellLists& split_cell_lists = split_cell_lists_levels_[level];
			for (size_t s = 0; s != split_cell_lists.size(); ++s) {
				StdLargeVec<CellList*>& cell_lists = split_cell_lists[s];
				parallel_for(blocked_range<size_t>(0, cell_lists.size()),
					[&](const blocked_range<size_t>& r) {
						for (size_t num = r.begin(); num < r.end(); ++num) {

							CellList* cell_list = cell_lists[num];
							int i = (int)cell_list->cell_location_[0];
							int j = (int)cell_list->cell_location_[1];
							ConcurrentListDataVector& particle_data_lists = cell_list->particle_data_lists_;
							for (size_t num = 0; num != cell_list->real_particle_count_; ++num) {

								ListData& list_data = particle_data_lists[num];
								size_t particle_index_here = list_data.first;
								Neighborhood& neighborhood_here = inner_configuration[particle_index_here];
								NeighborList& neighbor_list_here = std::get<0>(neighborhood_here);
								size_t previous_count_of_neigbors = std::get<2>(neighborhood_here);
								Vecu number_of_cells = number_of_cells_levels_[level];
								for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++l)
									for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++m) {
										ConcurrentListDataVector list_data_targets 
											= cell_linked_lists[l][m].particle_data_lists_;
										for (size_t n = 0; n != list_data_targets.size(); ++n)
										{
											size_t particle_index_there = list_data_targets[n].first;
											//displacement pointing from neighboring particle to origin particle
											Vecd displacement = list_data.second - list_data_targets[n].second;
											if (displacement.norm() < cell_spacing) {
												//neigbor particles for the original particle
												size_t current_count_here = std::get<1>(neighborhood_here);
												current_count_here >= previous_count_of_neigbors ?
													neighbor_list_here.push_back(
														new NeighborRelationWithVariableSmoothingLength(base_particle_data, *kernel_,
															displacement, particle_index_here, particle_index_there))
													: neighbor_list_here[current_count_here]->resetRelation(base_particle_data, *kernel_,
														displacement, particle_index_here, particle_index_there);
												std::get<1>(neighborhood_here)++;
												//neighbor particles for the target particle
												if (particle_index_there < base_particles_->real_particles_bound_) {
													Neighborhood& neighborhood_there = inner_configuration[particle_index_there];
													NeighborList& neighbor_list_there = std::get<0>(neighborhood_there);
													size_t current_count_there = std::get<1>(neighborhood_there);
													current_count_there >= std::get<2>(neighborhood_there) ?
														neighbor_list_there.push_back(
															new NeighborRelationWithVariableSmoothingLength(
																neighbor_list_here[current_count_here], particle_index_here))
														: neighbor_list_there[current_count_there]
														->resetSymmetricRelation(neighbor_list_here[current_count_here], particle_index_here);
													std::get<1>(neighborhood_there)++;
												}
											}
										}
									}
							}

						}
					}, ap);
			}
		}
	}
	//=================================================================================================//
	void MeshCellLinkedList
		::InsertACellLinkedListEntry(size_t particle_index, Vecd particle_position)
	{
		Vecu cellpos = GridIndexesFromPosition(particle_position);
		cell_linked_lists_[cellpos[0]][cellpos[1]].particle_data_lists_
			.push_back(make_pair(particle_index, particle_position));
	}
//=================================================================================================//
}
