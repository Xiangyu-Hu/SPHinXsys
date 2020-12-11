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
		BaseParticles* base_particles = sph_body_->base_particles_;
		Vecu number_of_cells = mesh_cell_linked_list_->NumberOfCells();
		matrix_cell cell_linked_lists = mesh_cell_linked_list_->CellLinkedLists();
		Kernel* current_kernel = sph_body_->kernel_;
		Real cutoff_radius_sqr = powern(current_kernel->GetCutOffRadius(), 2);

		parallel_for(blocked_range<size_t>(0, sph_body_->number_of_particles_),
			[&](const blocked_range<size_t>& r) {
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					Vecd particle_position = base_particles->pos_n_[num];
					Vecu cell_location = 
						mesh_cell_linked_list_->GridIndexFromPosition(particle_position);
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
									Vecd displacement = particle_position - target_particles[n].second;
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
}
