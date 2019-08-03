#include "mesh_cell_linked_list.h"
#include "base_body.h"
#include "sph_system.h"
#include "base_kernel.h"
#include "base_particles.h"


namespace SPH {

	MeshCellLinkedList::MeshCellLinkedList(Vecd lower_bound,
		Vecd upper_bound, Real cell_spacing, size_t buffer_size)
		: Mesh(lower_bound, upper_bound, cell_spacing, buffer_size)
	{
		number_of_cells_ = number_of_grid_points_;
		cutoff_radius_ = grid_spacing_;
		cell_spacing_ = grid_spacing_;
	}
	//===========================================================//
	int MeshCellLinkedList::ComputingSearchRage(int orign_refinement_level,
		int target_refinement_level)
	{
		return orign_refinement_level >= target_refinement_level
			? 1 : powern(2, target_refinement_level - orign_refinement_level);
	}
	//===========================================================//
	Kernel&  MeshCellLinkedList
		::ChoosingKernel(Kernel *orignal_kernel, Kernel *target_kernel)
	{
		return orignal_kernel->GetSmoothingLength() 
				> target_kernel->GetSmoothingLength()
			? *orignal_kernel : *target_kernel;
	}
	//===========================================================//
	void MeshCellLinkedList::UpdateParticleCellLocation(SPHBody &body)
	{
		StdVec<BaseParticleData> &base_particle_data = body.base_particles_->base_particle_data_;

		//rebuild the corresponding particle list.
		parallel_for(blocked_range<size_t>(0, base_particle_data.size()),
			[&](const blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				Vecu cellpos = GridIndexesFromPosition(base_particle_data[i].pos_n_);
				base_particle_data[i].cell_location_ = cellpos;
			}
		}, ap);
	}
	//===========================================================//
	void MeshCellLinkedList::ClearByCellParticleLists(SPHBody &body)
	{
		ByCellLists by_cell_lists_particle_indexes
			= body.by_cell_lists_particle_indexes_;

		for (size_t i = 0; i < body.number_of_by_cell_lists_; i++) {
			parallel_for(blocked_range<size_t>(0, by_cell_lists_particle_indexes[i].size()),
				[&](const blocked_range<size_t>& r) {
				for (size_t j = r.begin(); j != r.end(); ++j) 
					by_cell_lists_particle_indexes[i][j].clear();
			}, ap);
			by_cell_lists_particle_indexes[i].clear();
		}
	}
}