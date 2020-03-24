#include "mesh_cell_linked_list.h"
#include "array_allocation.h"
#include "sph_system.h"
#include "base_particles.h"
#include "base_body.h"
#include "base_kernel.h"
#include "neighboring_particle.h"
#include "base_data_package.h"

namespace SPH {

	CellList::CellList() : cell_location_(0)
	{
		particle_data_lists_.reserve(36);
	}
	//===============================================================//
	void MeshCellLinkedList
		::AllocateMeshDataMatrix()
	{
		Allocate3dArray(cell_linked_lists_, number_of_grid_points_);
		for (size_t i = 0; i != number_of_cells_[0]; ++i)
			for (size_t j = 0; j != number_of_cells_[1]; ++j)
				for (size_t k = 0; k != number_of_cells_[2]; ++k) {
					cell_linked_lists_[i][j][k].setCellInformation(Vecu(i, j, k));
				}
	}
	//===============================================================//
	void MeshCellLinkedList
		::DeleteMeshDataMatrix()
	{
		Delete3dArray(cell_linked_lists_, number_of_grid_points_);
	}
	//===============================================================//
	void MeshCellLinkedList::ClearCellLists()
	{
		parallel_for(blocked_range3d<size_t>(0, number_of_cells_[0], 0, number_of_cells_[1], 0, number_of_cells_[2]),
			[&](const blocked_range3d<size_t>& r) {
			for (size_t i = r.pages().begin(); i != r.pages().end(); ++i)
				for (size_t j = r.rows().begin(); j != r.rows().end(); ++j)
					for (size_t k = r.cols().begin(); k != r.cols().end(); ++k)
					{
						cell_linked_lists_[i][j][k].particle_data_lists_.clear();
						cell_linked_lists_[i][j][k].real_particles_in_cell_ = 0;

					}
		}, ap);
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateSplitCellLists(SPHBody &body)
	{
		//clear the data
		ClearSplitCellLists(body);

		SplitCellLists& split_cell_lists = body.split_cell_lists_;

		parallel_for(blocked_range<size_t>(0, split_cell_lists.size()),
			[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num) {
					/** convert 1d vector index to mesh index. */
					Vec3u split_mesh_index = transfer1DtoMeshIndex(Vec3u(3, 3, 3), num);
					size_t l = split_mesh_index[0];
					size_t m = split_mesh_index[1];
					size_t n = split_mesh_index[2];

					Vecu starting_index(2 - l, 2 - m, 2 - n);
					Vecu number_of_operation = (number_of_cells_ + starting_index) / 3;
					for (size_t i = 0; i != number_of_operation[0]; ++i)
						for (size_t j = 0; j != number_of_operation[1]; ++j)
							for (size_t k = 0; k != number_of_operation[2]; ++k) {
								CellList& cell_list = cell_linked_lists_[3 * i + l][3 * j + m][3 * k + n];
								size_t real_particles_in_cell = cell_list.particle_data_lists_.size();
								if (real_particles_in_cell != 0) {
									cell_list.real_particles_in_cell_ = real_particles_in_cell;
									split_cell_lists[num].push_back(&cell_linked_lists_[3 * i + l][3 * j + m][3 * k + n]);
								}
							}
				}
			}, ap);
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateInnerConfiguration(SPHBody &body,
		InnerParticleConfiguration& inner_particle_configuration)
	{
		StdLargeVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;

		parallel_for(blocked_range<size_t>(0, body.number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t num = r.begin(); num != r.end(); ++num)
			{
				Vecu cell_location = GridIndexesFromPosition(base_particle_data[num].pos_n_);
				int i = (int)cell_location[0];
				int j = (int)cell_location[1];
				int k = (int)cell_location[2];

				Neighborhood& neighborhood = inner_particle_configuration[num];
				NeighborList& neighbor_list = std::get<0>(neighborhood);
				size_t previous_count_of_neigbors = std::get<2>(neighborhood);

				for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
				{
					for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
					{
						for (int q = SMAX(k - 1, 0); q <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++q)
						{
							ConcurrentListDataVector &target_particles = cell_linked_lists_[l][m][q].particle_data_lists_;
							for (size_t n = 0; n != target_particles.size(); ++n)
							{
								//displacement pointing from neighboring particle to origin particle
								Vecd displacement = base_particle_data[num].pos_n_ - target_particles[n].second;
								if (displacement.norm() <= cutoff_radius_ && num != target_particles[n].first)
								{
									std::get<1>(neighborhood) >= previous_count_of_neigbors ?
										neighbor_list.push_back(new NeighboringParticle(*(body.kernel_), displacement, target_particles[n].first))
										: neighbor_list[std::get<1>(neighborhood)]->Reset(*(body.kernel_), displacement, target_particles[n].first);
									std::get<1>(neighborhood)++;
								}
							}
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
	//===============================================================//
	void MeshCellLinkedList::UpdateContactConfiguration(SPHBody &body)
	{
		StdLargeVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		ContatcParticleConfiguration& current_contact_configuration = body.current_contact_configuration_;
		ContactParticles& indexes_contact_particles = body.indexes_contact_particles_;
		SPHBodyVector contact_bodies = body.contact_map_.second;

		//clear previous surface particles
		for (size_t body_num = 0; body_num < contact_bodies.size(); ++body_num)
		{
			indexes_contact_particles[body_num].clear();
			MeshCellLinkedList &target_mesh_cell_linked_list
				= *(contact_bodies[body_num]->mesh_cell_linked_list_);
			Vecu target_number_of_cells = target_mesh_cell_linked_list.number_of_cells_;
			int search_range
				= ComputingSearchRage(body.refinement_level_,
					contact_bodies[body_num]->refinement_level_);
			Kernel &kernel = ChoosingKernel(body.kernel_, contact_bodies[body_num]->kernel_);
			Real cutoff_radius = kernel.GetCutOffRadius();

			parallel_for(blocked_range<size_t>(0, body.number_of_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					Vecu target_cell_index
						= target_mesh_cell_linked_list
						.GridIndexesFromPosition(base_particle_data[num].pos_n_);
					int i = (int)target_cell_index[0];
					int j = (int)target_cell_index[1];
					int k = (int)target_cell_index[2];

					matrix_cell target_cell_linked_lists
						= target_mesh_cell_linked_list.cell_linked_lists_;

					Neighborhood& neighborhood = current_contact_configuration[body_num][num];
					NeighborList& neighbor_list = std::get<0>(neighborhood);
					size_t previous_count_of_neigbors = std::get<2>(neighborhood);

					for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
							for (int q = SMAX(k - search_range, 0); q <= SMIN(k + search_range, int(target_number_of_cells[2]) - 1); ++q)
							{
								ConcurrentListDataVector &target_particles
									= target_cell_linked_lists[l][m][q].particle_data_lists_;
								for (size_t n = 0; n < target_particles.size(); n++)
								{
									//displacement pointing from neighboring particle to origin particle
									Vecd displacement = base_particle_data[num].pos_n_
										- target_particles[n].second;
									if (displacement.norm() <= cutoff_radius)
									{
										std::get<1>(neighborhood) >= previous_count_of_neigbors ?
											neighbor_list.push_back(new NeighboringParticle(kernel, displacement, target_particles[n].first))
											: neighbor_list[std::get<1>(neighborhood)]->Reset(kernel, displacement, target_particles[n].first);
										std::get<1>(neighborhood)++;
									}
								}
							}
					size_t current_count_of_neighbors = std::get<1>(neighborhood);
					neighbor_list.resize(current_count_of_neighbors);
					std::get<2>(neighborhood) = current_count_of_neighbors;
					std::get<1>(neighborhood) = 0;
					if (current_count_of_neighbors != 0)
						indexes_contact_particles[body_num].push_back(num);
				}
			}, ap);
		}
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateInteractionConfiguration(SPHBody &body,
		SPHBodyVector interacting_bodies)
	{
		StdLargeVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		ContatcParticleConfiguration& current_contact_configuration = body.current_contact_configuration_;
		ContactParticles& indexes_contact_particles = body.indexes_contact_particles_;
		SPHBodyVector contact_bodies = body.contact_map_.second;
		IndexVector contact_configuration_index(interacting_bodies.size());

		//clear previous interacting particles
		//build up contact configuration indexes
		for (size_t intertaction_body_num = 0;
			intertaction_body_num < interacting_bodies.size(); ++intertaction_body_num) {
			for (size_t contact_body_num = 0;
				contact_body_num < contact_bodies.size(); ++contact_body_num) {
				if (interacting_bodies[intertaction_body_num] == contact_bodies[contact_body_num]) {
					contact_configuration_index[intertaction_body_num] = contact_body_num;
					indexes_contact_particles[contact_body_num].clear();
				}
			}

			MeshCellLinkedList &target_mesh_cell_linked_list
				= *(interacting_bodies[intertaction_body_num]->mesh_cell_linked_list_);
			Vecu target_number_of_cells = target_mesh_cell_linked_list.number_of_cells_;
			int search_range
				= ComputingSearchRage(body.refinement_level_,
					interacting_bodies[intertaction_body_num]->refinement_level_);
			Kernel &kernel = ChoosingKernel(body.kernel_,
				interacting_bodies[intertaction_body_num]->kernel_);
			Real cutoff_radius = kernel.GetCutOffRadius();

			parallel_for(blocked_range<size_t>(0, body.number_of_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num) {
					Vecu target_cell_index
						= target_mesh_cell_linked_list
						.GridIndexesFromPosition(base_particle_data[num].pos_n_);
					int i = (int)target_cell_index[0];
					int j = (int)target_cell_index[1];
					int k = (int)target_cell_index[2];

					matrix_cell target_cell_linked_lists
						= interacting_bodies[intertaction_body_num]
						->mesh_cell_linked_list_->cell_linked_lists_;
					size_t contact_body_num
						= contact_configuration_index[intertaction_body_num];

					Neighborhood& neighborhood = current_contact_configuration[contact_body_num][num];
					NeighborList& neighbor_list = std::get<0>(neighborhood);
					size_t previous_count_of_neigbors = std::get<2>(neighborhood);

					for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
							for (int q = SMAX(k - search_range, 0); q <= SMIN(k + search_range, int(target_number_of_cells[2]) - 1); ++q)
							{
								ConcurrentListDataVector &target_particles
									= target_cell_linked_lists[l][m][q].particle_data_lists_;
								for (size_t n = 0; n < target_particles.size(); n++)
								{
									//displacement pointing from neighboring particle to origin particle
									Vecd displacement = base_particle_data[num].pos_n_
										- target_particles[n].second;
									if (displacement.norm() <= cutoff_radius)
									{
										std::get<1>(neighborhood) >= previous_count_of_neigbors ?
											neighbor_list.push_back(new NeighboringParticle(kernel, displacement, target_particles[n].first))
											: neighbor_list[std::get<1>(neighborhood)]->Reset(kernel, displacement, target_particles[n].first);
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
	//===============================================================//
	bool MeshCellLinkedList::IsInSequence(Vecd pos_1, Vecd pos_2)
	{
		size_t index_1 = clamp((int)floor(pos_1[0] / cell_spacing_), 0, int(number_of_cells_[0]) - 1)
			*number_of_cells_[1] * number_of_cells_[2]
			+ clamp((int)floor(pos_1[1] / cell_spacing_), 0, int(number_of_cells_[1]) - 1)*number_of_cells_[2]
			+ clamp((int)floor(pos_1[2] / cell_spacing_), 0, int(number_of_cells_[2]) - 1);
		size_t index_2 = clamp((int)floor(pos_2[0] / cell_spacing_), 0, int(number_of_cells_[0]) - 1)
			*number_of_cells_[1] * number_of_cells_[2]
			+ clamp((int)floor(pos_2[1] / cell_spacing_), 0, int(number_of_cells_[1]) - 1)*number_of_cells_[2]
			+ clamp((int)floor(pos_2[2] / cell_spacing_), 0, int(number_of_cells_[2]) - 1);

		return index_1 > index_2 ? false : true;
	}
	//===========================================================//
	void MeshCellLinkedList
		::InsertACellLinkedListEntry(size_t particle_index, Vecd particle_position, Vecu cellpos)
	{
		cell_linked_lists_[cellpos[0]][cellpos[1]][cellpos[2]].particle_data_lists_
			.push_back(make_pair(particle_index, particle_position));
	}
	//===========================================================//
}
