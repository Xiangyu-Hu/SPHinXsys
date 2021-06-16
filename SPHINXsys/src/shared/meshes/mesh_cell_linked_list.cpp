/**
 * @file 	mesh_cell_linked_list.cpp
 * @author	Yongchuan Yu, Chi ZHang and Xiangyu Hu
 */

#include "mesh_cell_linked_list.h"
#include "base_kernel.h"
#include "base_body.h"
#include "particle_adaptation.h"
#include "base_particles.h"


namespace SPH {
	//=================================================================================================//
	BaseMeshCellLinkedList::
		BaseMeshCellLinkedList(SPHBody& sph_body, ParticleAdaptation& particle_adaptation,
			BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width) : 
		Mesh(tentative_bounds, grid_spacing, buffer_width),
		sph_body_(sph_body), kernel_(*particle_adaptation.getKernel()), 
		base_particles_(NULL) {}
	//=================================================================================================//
	void BaseMeshCellLinkedList::clearSplitCellLists(SplitCellLists& split_cell_lists)
	{
		for (size_t i = 0; i < split_cell_lists.size(); i++)
			split_cell_lists[i].clear();
	}
	//=================================================================================================//
	MeshCellLinkedList::MeshCellLinkedList(SPHBody& sph_body, ParticleAdaptation& particle_adaptation, 
		BoundingBox tentative_bounds, Real grid_spacing, size_t buffer_width)
		: BaseMeshCellLinkedList(sph_body, particle_adaptation, tentative_bounds, grid_spacing, buffer_width)
	{
		name_ = "MeshCellLinkedList";
		allocateMeshDataMatrix();
	}
	//=================================================================================================//
	void MeshCellLinkedList::UpdateCellLists()
	{
		clearCellLists();
		StdLargeVec<Vecd>& pos_n = base_particles_->pos_n_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		parallel_for(blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					insertACellLinkedParticleIndex(i, pos_n[i]);
				}
			}, ap);
		UpdateCellListData();
		updateSplitCellLists(sph_body_.split_cell_lists_);
	}
	//=================================================================================================//
	void MeshCellLinkedList::assignBaseParticles(BaseParticles* base_particles)
	{
		base_particles_ = base_particles;
	};
	//=================================================================================================//
	void MeshCellLinkedList::computingSequence(StdLargeVec<size_t>& sequence)
	{
		StdLargeVec<Vecd>& positions = base_particles_->pos_n_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		parallel_for(blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
                    sequence[i] = transferMeshIndexToMortonOrder(CellIndexFromPosition(positions[i]));
				}
			}, ap);
	}
	//=================================================================================================//
	MultilevelMeshCellLinkedList
		::MultilevelMeshCellLinkedList(SPHBody& sph_body, ParticleAdaptation& particle_adaptation, 
			BoundingBox tentative_bounds, Real reference_grid_spacing, 
			size_t total_levels, Real maximum_spacing_ratio) :
		MultilevelMesh<SPHBody, BaseMeshCellLinkedList, MeshCellLinkedList>(sph_body, particle_adaptation, 
			tentative_bounds, reference_grid_spacing, total_levels, maximum_spacing_ratio, 2),
			h_ratio_(dynamic_cast<ParticleWithLocalRefinement&>(particle_adaptation).h_ratio_)
	{
		name_ = "MultilevelMeshCellLinkedList";
	}
	//=================================================================================================//
	size_t MultilevelMeshCellLinkedList::getMeshLevel(Real particle_cutoff_radius)
	{
		for (size_t level = total_levels_; level != 0; --level)
			if (particle_cutoff_radius - mesh_levels_[level - 1]->GridSpacing() < Eps) return level - 1; //jump out the loop!

		std::cout << "\n Error: MeshCellLinkedList level searching out of bound!" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
		return 999; //means an error in level searching
	};
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::
		insertACellLinkedParticleIndex(size_t particle_index, Vecd particle_position)
	{
		size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
		mesh_levels_[level]->insertACellLinkedParticleIndex(particle_index, particle_position);
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::
		InsertACellLinkedListDataEntry(size_t particle_index, Vecd particle_position)
	{
		size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
		mesh_levels_[level]->InsertACellLinkedListDataEntry(particle_index, particle_position);
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::assignBaseParticles(BaseParticles* base_particles)
	{
		base_particles_ = base_particles;
		for (size_t l = 0; l != total_levels_; ++l) {
			mesh_levels_[l]->assignBaseParticles(base_particles);
		}
	};
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::UpdateCellLists()
	{
		for (size_t level = 0; level != total_levels_; ++level) 
			mesh_levels_[level]->clearCellLists();
	
		StdLargeVec<Vecd>& pos_n = base_particles_->pos_n_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		//rebuild the corresponding particle list.
		parallel_for(blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i) {
					insertACellLinkedParticleIndex(i, pos_n[i]);
				}
			}, ap);

		for (size_t level = 0; level != total_levels_; ++level) 
			mesh_levels_[level]->UpdateCellListData();
		updateSplitCellLists(sph_body_.split_cell_lists_);
	}
	//=================================================================================================//
	void MultilevelMeshCellLinkedList::
		tagBodyPartByCell(CellLists& cell_lists, std::function<bool(Vecd, Real)>& check_included)
	{
		for (size_t l = 0; l != total_levels_; ++l) {
			mesh_levels_[l]->tagBodyPartByCell(cell_lists, check_included);
		}
	}
	//=================================================================================================//
}
