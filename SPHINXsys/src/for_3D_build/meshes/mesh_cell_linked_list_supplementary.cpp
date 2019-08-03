#include "mesh_cell_linked_list.h"
#include "array_allocation.h"
#include "sph_system.h"
#include "base_particles.h"
#include "base_body.h"
#include "base_kernel.h"
#include "neighboring_particle.h"
#include "base_data_package.h"

namespace SPH {

	CellList::CellList()
	{
		particle_data_lists_.reserve(27);

	}
	//===============================================================//
	void MeshCellLinkedList
		::AllocateMeshDataMatrix()
	{
		Allocate3dArray(cell_linked_lists_, number_of_grid_points_);
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
					}
		}, ap);
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateCellLists(SPHBody &body)
	{
		ClearCellLists();
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;

		//rebuild the corresponding particle list.
		parallel_for(blocked_range<size_t>(0, base_particle_data.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i != r.end(); ++i)
			{
				Vecu cellpos = GridIndexesFromPosition(base_particle_data[i].pos_n_);
				base_particle_data[i].cell_location_ = cellpos;
				cell_linked_lists_[cellpos[0]][cellpos[1]][cellpos[2]]
					.particle_data_lists_.push_back(make_pair(i, base_particle_data[i].pos_n_));
			}
		}, ap);
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateByCellParticleLists(SPHBody &body)
	{
		cout << "\n The function " 
			<< "MeshCellLinkedList::UpdateByCellParticleLists"
			<< " is not done in 3D. Exit the program! \n";
		exit(0);
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateInnerConfiguration(SPHBody &body)
	{
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		NeighborList current_inner_configuration = body.current_inner_configuration_;

		parallel_for(blocked_range<size_t>(0, body.number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t num = r.begin(); num != r.end(); ++num)
			{
				int i = (int)base_particle_data[num].cell_location_[0];
				int j = (int)base_particle_data[num].cell_location_[1];
				int k = (int)base_particle_data[num].cell_location_[2];

				StdVec<NeighboringParticle> &current_neighbor_list = current_inner_configuration[num];
				size_t current_number_of_neighbors = current_neighbor_list.size();
				size_t count_of_neigbors = 0;

				for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
				{
					for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
					{
						for (int q = SMAX(k - 1, 0); q <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++q)
						{
							ListDataVector &target_particles = cell_linked_lists_[l][m][q].particle_data_lists_;
							for (size_t n = 0; n != target_particles.size(); ++n)
							{
								//displacement pointing from neighboring particle to origin particle
								Vecd displacement = base_particle_data[num].pos_n_ - target_particles[n].second;
								if (displacement.norm() <= cutoff_radius_ && num != target_particles[n].first)
								{
									count_of_neigbors >= current_number_of_neighbors ?
										current_neighbor_list
										.push_back(NeighboringParticle(*(body.kernel_), displacement, target_particles[n].first))
										: current_neighbor_list[count_of_neigbors]
										.Reset(*(body.kernel_), displacement, target_particles[n].first);
									count_of_neigbors++;
								}
							}
						}
					}
				}
				//std::cout<<"count_of_neigbors :"<<count_of_neigbors<<std::endl;
				current_neighbor_list.resize(count_of_neigbors);
			}
		}, ap);
	}
	//===============================================================//
	void MeshCellLinkedList
		::BuildReferenceInnerConfiguration(SPHBody &body)
	{
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		ReferenceNeighborList reference_inner_configuration = body.reference_inner_configuration_;

		parallel_for(blocked_range<size_t>(0, body.number_of_particles_),
			[&](const blocked_range<size_t>& r) {
			for (size_t num = r.begin(); num != r.end(); ++num)
			{
				int i = (int)base_particle_data[num].cell_location_[0];
				int j = (int)base_particle_data[num].cell_location_[1];
				int k = (int)base_particle_data[num].cell_location_[2];

				StdVec<ReferenceNeighboringParticle> &reference_neighbor_list = reference_inner_configuration[num];
				size_t current_number_of_neighbors = reference_neighbor_list.size();
				size_t count_of_neigbors = 0;

				for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells_[0]) - 1); ++l)
				{
					for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells_[1]) - 1); ++m)
					{
						for (int q = SMAX(k - 1, 0); q <= SMIN(k + 1, int(number_of_cells_[2]) - 1); ++q)
						{
							ListDataVector &target_particles = cell_linked_lists_[l][m][q].particle_data_lists_;
							for (size_t n = 0; n != target_particles.size(); ++n)
							{
								//displacement pointing from neighboring particle to origin particle
								Vecd displacement = base_particle_data[num].pos_n_ - target_particles[n].second;
								if (displacement.norm() <= cutoff_radius_ && num != target_particles[n].first)
								{
									count_of_neigbors >= current_number_of_neighbors ?
										reference_neighbor_list
										.push_back(ReferenceNeighboringParticle(*(body.kernel_), displacement, target_particles[n].first))
										: reference_neighbor_list[count_of_neigbors]
										.Reset(*(body.kernel_), displacement, target_particles[n].first);
									count_of_neigbors++;
								}
							}
						}
					}
				}
				//std::cout<<"count_of_neigbors :"<<count_of_neigbors<<std::endl;
				reference_neighbor_list.resize(count_of_neigbors);
			}
		}, ap);
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateContactConfiguration(SPHBody &body)
	{
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		ContactNeighborList  current_contact_configuration = body.current_contact_configuration_;
		ContactParticleList indexes_contact_particles = body.indexes_contact_particles_;
		StdVec<SPHBody*> contact_bodies = body.contact_map_.second;

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
					StdVec<NeighboringParticle> &current_neighbor_list
						= current_contact_configuration[body_num][num];
					size_t current_number_of_neighbors = current_neighbor_list.size();
					size_t count_of_neigbors = 0;

					for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
							for (int q = SMAX(k - search_range, 0); q <= SMIN(k + search_range, int(target_number_of_cells[2]) - 1); ++q)
							{
								ListDataVector &target_particles
									= target_cell_linked_lists[l][m][q].particle_data_lists_;
								for (size_t n = 0; n < target_particles.size(); n++)
								{
									//displacement pointing from neighboring particle to origin particle
									Vecd displacement = base_particle_data[num].pos_n_
										- target_particles[n].second;
									if (displacement.norm() <= cutoff_radius)
									{
										count_of_neigbors >= current_number_of_neighbors ?
											current_neighbor_list
											.push_back(NeighboringParticle(kernel, displacement, target_particles[n].first))
											: current_neighbor_list[count_of_neigbors]
											.Reset(kernel, displacement, target_particles[n].first);
										count_of_neigbors++;
									}
								}
							}
					current_neighbor_list.resize(count_of_neigbors);
					if (current_neighbor_list.size() != 0)
					{
						indexes_contact_particles[body_num].push_back(num);
					}
				}
			}, ap);
		}
	}
	//===============================================================//
	void MeshCellLinkedList
		::BuildReferenceContactConfiguration(SPHBody &body)
	{
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		ReferenceContactNeighborList  reference_contact_configuration = body.reference_contact_configuration_;
		ContactParticleList indexes_contact_particles = body.indexes_contact_particles_;
		StdVec<SPHBody*> contact_bodies = body.contact_map_.second;

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
					StdVec<ReferenceNeighboringParticle> &reference_neighbor_list
						= reference_contact_configuration[body_num][num];
					size_t current_number_of_neighbors = reference_neighbor_list.size();
					size_t count_of_neigbors = 0;

					for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
							for (int q = SMAX(k - search_range, 0); q <= SMIN(k + search_range, int(target_number_of_cells[2]) - 1); ++q)
							{
								ListDataVector &target_particles
									= target_cell_linked_lists[l][m][q].particle_data_lists_;
								for (size_t n = 0; n < target_particles.size(); n++)
								{
									//displacement pointing from neighboring particle to origin particle
									Vecd displacement = base_particle_data[num].pos_n_
										- target_particles[n].second;
									if (displacement.norm() <= cutoff_radius)
									{
										count_of_neigbors >= current_number_of_neighbors ?
											reference_neighbor_list
											.push_back(ReferenceNeighboringParticle(kernel, displacement, target_particles[n].first))
											: reference_neighbor_list[count_of_neigbors]
											.Reset(kernel, displacement, target_particles[n].first);
										count_of_neigbors++;
									}
								}
							}
					reference_neighbor_list.resize(count_of_neigbors);
					if (reference_neighbor_list.size() != 0)
					{
						indexes_contact_particles[body_num].push_back(num);
					}
				}
			}, ap);
		}
	}
	//===============================================================//
	void MeshCellLinkedList::UpdateInteractionConfiguration(SPHBody &body,
		SPHBodyVector interacting_bodies)
	{
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;
		ContactNeighborList  current_contact_configuration = body.current_contact_configuration_;
		ContactParticleList indexes_contact_particles = body.indexes_contact_particles_;
		StdVec<SPHBody*> contact_bodies = body.contact_map_.second;
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
					StdVec<NeighboringParticle> &current_neighbor_list
						= current_contact_configuration[contact_body_num][num];
					size_t current_number_of_neighbors = current_neighbor_list.size();
					size_t count_of_neigbors = 0;

					for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
						for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
							for (int q = SMAX(k - search_range, 0); q <= SMIN(k + search_range, int(target_number_of_cells[2]) - 1); ++q)
							{
								ListDataVector &target_particles
									= target_cell_linked_lists[l][m][k].particle_data_lists_;
								for (size_t n = 0; n < target_particles.size(); n++)
								{
									//displacement pointing from neighboring particle to origin particle
									Vecd displacement = base_particle_data[num].pos_n_
										- target_particles[n].second;
									if (displacement.norm() <= cutoff_radius)
									{
										count_of_neigbors >= current_number_of_neighbors ?
											current_neighbor_list
											.push_back(NeighboringParticle(kernel, displacement, target_particles[n].first))
											: current_neighbor_list[count_of_neigbors]
											.Reset(kernel, displacement, target_particles[n].first);
										count_of_neigbors++;
									}
								}
							}
					current_neighbor_list.resize(count_of_neigbors);
					if (current_neighbor_list.size() != 0)
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
}
