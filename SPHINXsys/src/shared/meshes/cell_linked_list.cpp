/**
 * @file 	cell_linked_list.cpp
 * @author	Yongchuan Yu, Chi ZHang and Xiangyu Hu
 */

#include "cell_linked_list.h"
#include "base_kernel.h"
#include "base_body.h"
#include "adaptation.h"
#include "base_particles.h"
#include "base_particle_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	BaseCellLinkedList::
		BaseCellLinkedList(RealBody &real_body, SPHAdaptation &sph_adaptation)
		: BaseMeshField("CellLinkedList"),
		  real_body_(real_body), kernel_(*sph_adaptation.getKernel()),
		  base_particles_(nullptr) {}
	//=================================================================================================//
	void BaseCellLinkedList::clearSplitCellLists(SplitCellLists &split_cell_lists)
	{
		for (size_t i = 0; i < split_cell_lists.size(); i++)
			split_cell_lists[i].clear();
	}
	//=================================================================================================//
	CellLinkedList::CellLinkedList(BoundingBox tentative_bounds, Real grid_spacing,
								   RealBody &real_body, SPHAdaptation &sph_adaptation)
		: BaseCellLinkedList(real_body, sph_adaptation), Mesh(tentative_bounds, grid_spacing, 2)
	{
		allocateMeshDataMatrix();
	}
	//=================================================================================================//
	void CellLinkedList::UpdateCellLists()
	{
		clearCellLists();
		StdLargeVec<Vecd> &pos_n = base_particles_->pos_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		parallel_for(
			blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					insertParticleIndex(i, pos_n[i]);
				}
			},
			ap);

		UpdateCellListData();

		if (real_body_.getUseSplitCellLists())
		{
			updateSplitCellLists(real_body_.getSplitCellLists());
		}
	}
	//=================================================================================================//
	void CellLinkedList::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
	};
	//=================================================================================================//
	void CellLinkedList::computingSequence(StdLargeVec<size_t> &sequence)
	{
		StdLargeVec<Vecd> &positions = base_particles_->pos_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		parallel_for(
			blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					sequence[i] = transferMeshIndexToMortonOrder(CellIndexFromPosition(positions[i]));
				}
			},
			ap);
	}
	//=================================================================================================//
	MultilevelCellLinkedList::
		MultilevelCellLinkedList(BoundingBox tentative_bounds, Real reference_grid_spacing,
								 size_t total_levels, RealBody &real_body, SPHAdaptation &sph_adaptation)
		: MultilevelMesh<BaseCellLinkedList, CellLinkedList, RefinedMesh<CellLinkedList>>(
			  tentative_bounds, reference_grid_spacing, total_levels, real_body, sph_adaptation),
		  h_ratio_(DynamicCast<ParticleWithLocalRefinement>(this, &sph_adaptation)->h_ratio_) {}
	//=================================================================================================//
	size_t MultilevelCellLinkedList::getMeshLevel(Real particle_cutoff_radius)
	{
		for (size_t level = total_levels_; level != 0; --level)
			if (particle_cutoff_radius - mesh_levels_[level - 1]->GridSpacing() < Eps)
				return level - 1; // jump out the loop!

		std::cout << "\n Error: CellLinkedList level searching out of bound!" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
		return 999; // means an error in level searching
	};
	//=================================================================================================//
	void MultilevelCellLinkedList::
		insertParticleIndex(size_t particle_index, const Vecd &particle_position)
	{
		size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
		mesh_levels_[level]->insertParticleIndex(particle_index, particle_position);
	}
	//=================================================================================================//
	void MultilevelCellLinkedList::
		InsertListDataEntry(size_t particle_index, const Vecd &particle_position)
	{
		size_t level = getMeshLevel(kernel_.CutOffRadius(h_ratio_[particle_index]));
		mesh_levels_[level]->InsertListDataEntry(particle_index, particle_position);
	}
	//=================================================================================================//
	void MultilevelCellLinkedList::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
		for (size_t l = 0; l != total_levels_; ++l)
		{
			mesh_levels_[l]->assignBaseParticles(base_particles);
		}
	};
	//=================================================================================================//
	void MultilevelCellLinkedList::UpdateCellLists()
	{
		for (size_t level = 0; level != total_levels_; ++level)
			mesh_levels_[level]->clearCellLists();

		StdLargeVec<Vecd> &pos_n = base_particles_->pos_;
		size_t total_real_particles = base_particles_->total_real_particles_;
		// rebuild the corresponding particle list.
		parallel_for(
			blocked_range<size_t>(0, total_real_particles),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					insertParticleIndex(i, pos_n[i]);
				}
			},
			ap);

		for (size_t level = 0; level != total_levels_; ++level)
		{
			mesh_levels_[level]->UpdateCellListData();
		}

		if (real_body_.getUseSplitCellLists())
		{
			updateSplitCellLists(real_body_.getSplitCellLists());
		}
	}
	//=================================================================================================//
	void MultilevelCellLinkedList::
		tagBodyPartByCell(ConcurrentIndexesInCells &cell_lists, std::function<bool(Vecd, Real)> &check_included)
	{
		for (size_t l = 0; l != total_levels_; ++l)
		{
			mesh_levels_[l]->tagBodyPartByCell(cell_lists, check_included);
		}
	}
	//=================================================================================================//
}
