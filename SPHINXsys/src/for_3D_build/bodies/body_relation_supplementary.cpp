/**
 * @file 	body_relation_supplementary.cpp
 * @brief 	Here, Functions belong to BaseBody, RealBody and FictitiousBody are given.
 * @author	hi ZHang and Xiangyu Hu
 * @version	0.1
 * 			0.2.0
 * 			Cell splitting algorithm are added.
 * 			Chi Zhang
 */
#include "body_relation.h"
#include "base_particles.h"
#include "base_kernel.h"
#include "mesh_cell_linked_list.h"
namespace SPH
{
	//=================================================================================================//
	void SPHBodyInnerRelation::updateConfiguration()
	{
		StdLargeVec<BaseParticleData>& base_particle_data
			= body_->base_particles_->base_particle_data_;
		Vecu number_of_cells = base_mesh_cell_linked_list_->NumberOfCells();
		matrix_cell cell_linked_lists = base_mesh_cell_linked_list_->CellLinkedLists();
		Kernel* current_kernel = body_->kernel_;
		Real cutoff_radius_sqr = powern(current_kernel->GetCutOffRadius(), 2);

		parallel_for(blocked_range<size_t>(0, body_->number_of_particles_),
			[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					BaseParticleData& base_particle_data_i = base_particle_data[num];
					Vecu cell_location = 
						base_mesh_cell_linked_list_->GridIndexFromPosition(base_particle_data_i.pos_n_);
					int i = (int)cell_location[0];
					int j = (int)cell_location[1];
					int k = (int)cell_location[2];

					Neighborhood& neighborhood = inner_configuration_[num];
					size_t current_count_of_neighbors = 0;
					for (int l = SMAX(i - 1, 0); l <= SMIN(i + 1, int(number_of_cells[0]) - 1); ++l)
					{
						for (int m = SMAX(j - 1, 0); m <= SMIN(j + 1, int(number_of_cells[1]) - 1); ++m)
						{
							for (int q = SMAX(k - 1, 0); q <= SMIN(k + 1, int(number_of_cells[2]) - 1); ++q)
							{
								CellListDataVector& target_particles = cell_linked_lists[l][m][q].cell_list_data_;
								for (size_t n = 0; n != target_particles.size(); ++n)
								{
									//displacement pointing from neighboring particle to origin particle
									Vecd displacement = base_particle_data_i.pos_n_ - target_particles[n].second;
									if (displacement.normSqr() <= cutoff_radius_sqr && num != target_particles[n].first)
									{
										current_count_of_neighbors >= neighborhood.memory_size_ ?
											createNeighborRelation(neighborhood,
												*current_kernel, displacement, num, target_particles[n].first)
											: initializeNeighborRelation(neighborhood, current_count_of_neighbors,
												*current_kernel, displacement, num, target_particles[n].first);
										current_count_of_neighbors++;
									}
								}
							}
						}
					}
					neighborhood.current_size_ = current_count_of_neighbors;
				}
			}, ap);
	}
	//=================================================================================================//
	void SPHBodyContactRelation::updateConfiguration()
	{
		StdLargeVec<BaseParticleData>& base_particle_data = body_->base_particles_->base_particle_data_;

		for (size_t relation_body_num = 0; relation_body_num < relation_bodies_.size(); ++relation_body_num) 
		{
			StdLargeVec<BaseParticleData>& target_base_particle_data
				= relation_bodies_[relation_body_num]->base_particles_->base_particle_data_;
			BaseMeshCellLinkedList& target_mesh_cell_linked_list
				= *(relation_bodies_[relation_body_num]->base_mesh_cell_linked_list_);
			Vecu target_number_of_cells = target_mesh_cell_linked_list.NumberOfCells();
			int search_range =
				base_mesh_cell_linked_list_->ComputingSearchRage(body_->refinement_level_,
					relation_bodies_[relation_body_num]->refinement_level_);
			Kernel& current_kernel = base_mesh_cell_linked_list_->ChoosingKernel(body_->kernel_,
				relation_bodies_[relation_body_num]->kernel_);
			Real cutoff_radius_sqr = powern(current_kernel.GetCutOffRadius(), 2);

			parallel_for(blocked_range<size_t>(0, body_->number_of_particles_),
				[&](const blocked_range<size_t>& r) {
					for (size_t num = r.begin(); num != r.end(); ++num) {

						BaseParticleData& base_particle_data_i = base_particle_data[num];
						Vecu target_cell_index = target_mesh_cell_linked_list
							.GridIndexFromPosition(base_particle_data_i.pos_n_);
						int i = (int)target_cell_index[0];
						int j = (int)target_cell_index[1];
						int k = (int)target_cell_index[2];

						matrix_cell target_cell_linked_lists
							= target_mesh_cell_linked_list.CellLinkedLists();

						Neighborhood& neighborhood = contact_configuration_[relation_body_num][num];
						size_t current_count_of_neighbors = 0;
						for (int l = SMAX(i - search_range, 0); l <= SMIN(i + search_range, int(target_number_of_cells[0]) - 1); ++l)
							for (int m = SMAX(j - search_range, 0); m <= SMIN(j + search_range, int(target_number_of_cells[1]) - 1); ++m)
								for (int q = SMAX(k - search_range, 0); q <= SMIN(k + search_range, int(target_number_of_cells[2]) - 1); ++q)
								{
									CellListDataVector& target_particles = target_cell_linked_lists[l][m][q].cell_list_data_;
									for (size_t n = 0; n < target_particles.size(); n++)
									{
										BaseParticleData& base_particle_data_j = target_base_particle_data[target_particles[n].first];
										//displacement pointing from neighboring particle to origin particle
										Vecd displacement = base_particle_data_i.pos_n_ - target_particles[n].second;
										if (checkNeighbor(displacement.normSqr(), cutoff_radius_sqr, base_particle_data_i, base_particle_data_j))
										{
											current_count_of_neighbors >= neighborhood.memory_size_ ?
												createNeighborRelation(neighborhood,
													current_kernel, displacement, num, target_particles[n].first)
												: initializeNeighborRelation(neighborhood, current_count_of_neighbors,
													current_kernel, displacement, num, target_particles[n].first);
											current_count_of_neighbors++;
										}
									}
								}
						neighborhood.current_size_ = current_count_of_neighbors;
					}
				}, ap);
		}
	}
	//=================================================================================================//
}
